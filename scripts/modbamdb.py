#!/usr/bin/env python

import os
import sys
import pysam
import argparse
import logging
import sqlite3

FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def rc(dna):
    ''' reverse complement '''
    complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def split_ml(mod_strings, ml):
    mls = []

    total_ms = sum([len(ms.split(',')[1:]) for ms in mod_strings])

    assert total_ms == len(ml), 'mod bam formatting error'

    i = 0
    for mod_string in mod_strings:
        m = mod_string.split(',')[1:] # discard first item (desc of mod base)
        mls.append(ml[i:i+len(m)])
        i += len(m)

        assert len(m) == len(mls[-1]), 'mod bam formatting error'
    
    return mls


def main(args):
    bam = pysam.AlignmentFile(args.bam)

    db_fn = '.'.join(args.bam.split('.')[:-1]) + '.db'

    if args.db:
        if args.db.endswith('.db'):
            db_fn = args.db
        else:
            db_fn = args.db + '.db'

    logger.info('db filename: %s' % db_fn)

    if os.path.exists(db_fn) and not args.append:
        sys.exit('database %s already exists' % db_fn)

    if args.append and not os.path.exists(db_fn):
        sys.exit('database %s does not exist and --append has been called' % db_fn)

    conn = sqlite3.connect(db_fn)

    c = conn.cursor()

    if not args.append:
        c.execute('''CREATE TABLE methdata (chrom text, pos integer, strand text, readname text, stat real, methstate integer, modname text)''')
        c.execute('''CREATE TABLE modnames (mod text)''')
        c.execute('''CREATE TABLE cutoffs (upper real, lower real, modname text)''')
        c.execute('''CREATE INDEX read_index ON methdata(readname)''')

    else:
        logger.info('appending records to %s' % db_fn)
    
    progress_interval = 1000000
    mod_count = 0
    ins_data = []
    minprob = float(args.minprob)
    mod_types = []

    for rec in bam.fetch(until_eof=True):
        if rec.is_unmapped:
            continue

        ap = dict([(k, v) for (k, v) in rec.get_aligned_pairs() if k is not None])

        mm = str(rec.get_tag('Mm')).rstrip(';')

        try:
            ml = rec.get_tag('Ml')
        except KeyError:
            continue

        mod_strings = mm.split(';')
        mls = split_ml(mod_strings, ml)

        seq = rec.seq

        read_str = '+'

        if rec.is_reverse:
            seq = rc(seq)
            read_str = '-'

        for mod_string, scores in zip(mod_strings, mls):
            m = mod_string.split(',')

            mod_info = m[0]

            mod_relpos = list(map(int, m[1:]))

            mod_strand = '+'

            if '-' in mod_info:
                mod_strand = '-'

            try:
                mod_base, mod_type = mod_info.split(mod_strand)
                mod_type = mod_type.rstrip('?.')
            except ValueError:
                logger.debug('invalid mod info for read %s' % rec.qname)
                continue

            assert len(mod_type) == 1, 'multiple modfications listed this way: %s is not yet supported, please send me an example!' % mod_info

            if mod_type not in mod_types:
                mod_types.append(mod_type)

            base_pos = [i for i, b in enumerate(seq) if b == mod_base]

            i = -1

            for skip, score in zip(mod_relpos, scores):
                mod_count += 1
                i += 1

                if skip > 0:
                    i += skip

                genome_pos = ap[base_pos[i]]

                if rec.is_reverse:
                    genome_pos = ap[len(rec.seq)-base_pos[i]-1]
                    # adjust position of - strand calls based on mod motif size (default = 2 as CG is probably the most frequent use case)
                    genome_pos -= (int(args.motifsize)-1)

                mod_prob = score/255
                can_prob = 1-mod_prob

                assert mod_prob <= 1.0

                methcall = 0

                if mod_prob >= minprob:
                    methcall = 1

                if can_prob >= minprob:
                    methcall = -1

                ins_data.append((rec.reference_name, genome_pos, read_str, rec.qname, mod_prob, methcall, mod_type))

                if mod_count % progress_interval == 0:
                    conn.executemany('INSERT INTO methdata VALUES (?,?,?,?,?,?,?)', ins_data)
                    ins_data = []
                    logger.info('processed %d records from %s' % (mod_count, args.bam))

        if len(ins_data) > 0:
            conn.executemany('INSERT INTO methdata VALUES (?,?,?,?,?,?,?)', ins_data)
            ins_data = []

    if len(ins_data) > 0:
        conn.executemany('INSERT INTO methdata VALUES (?,?,?,?,?,?,?)', ins_data)
        ins_data = []

    for mod in mod_types:
        c.execute("INSERT INTO modnames VALUES ('%s')" % mod)
        c.execute("INSERT INTO cutoffs VALUES ('%.4f', '%.4f', '%s')" % (minprob, 1-minprob, mod))

    logger.info('commiting %d records from %s to %s' % (mod_count, args.bam, db_fn))
    conn.commit()

    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='transfer modified base annotations from a methylartist .db to a .bam')
    parser.add_argument('-b', '--bam', required=True, help='bam used for methylation calling')
    parser.add_argument('-d', '--db', default=None, help='database name (default: auto-infer)')
    parser.add_argument('-p', '--minprob', default=0.8, help='probability threshold for calling modified or unmodified base (default = 0.8)')
    parser.add_argument('-a', '--append', default=False, action='store_true', help='append to database')
    parser.add_argument('--motifsize', default=2, help='mod motif size (default is 2 as "CG" is most common use case, e.g. set to 1 for 6mA)')
    args = parser.parse_args()
    main(args)