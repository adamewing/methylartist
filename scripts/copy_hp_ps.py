#!/usr/bin/env python

import argparse
import pysam
from collections import defaultdict as dd

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def main(args):
    assert args.haplobam.endswith('.bam')
    assert args.targetbam.endswith('.bam')

    out_fn = '.'.join(args.targetbam.split('.')[:-1]) + '.copytags.bam'

    if args.outbam:
        assert args.outbam.endswith('.bam')
        out_fn = args.outbam
    
    haplobam = pysam.AlignmentFile(args.haplobam)
    targetbam = pysam.AlignmentFile(args.targetbam)
    outbam = pysam.AlignmentFile(out_fn, 'wb', template=targetbam)

    data = dd(list)

    logger.info('collecting tagged reads...')

    tag_count = 0

    for i, read in enumerate(haplobam.fetch()):
        if read.has_tag('PS') and read.has_tag('HP'):
            data[read.qname].append([read.reference_name, read.pos, read.query_length, read.get_tag('PS'), read.get_tag('HP')])
            tag_count += 1
        
        if i % 100000 == 0:
            logger.info('parsed %d reads, %d tagged' % (i, tag_count))

    logger.info('collected %d tagged reads from %d reads' % (tag_count, i))

    tag_count = 0

    for i, read in enumerate(targetbam.fetch()):
        if read.qname in data:
            for rec in data[read.qname]:
                chrom, pos, qlen, ps, hp = rec
                if chrom != read.reference_name:
                    continue

                if read.pos < pos-qlen or read.pos > pos+qlen:
                    continue

                if read.has_tag('HP') or read.has_tag('PS'):
                    continue

                read.set_tags([("PS",ps),("HP",hp)])
                tag_count += 1

        outbam.write(read)

        if i % 100000 == 0:
            logger.info('wrote %d reads, %d tagged' % (i, tag_count))

    logger.info('wrote %d reads to %s, %d tagged' % (i, out_fn, tag_count))

    outbam.close()
    targetbam.close()
    haplobam.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='transfer modified base annotations from a methylartist .db to a .bam')
    parser.add_argument('--haplobam', required=True, help='haplotagged bam')
    parser.add_argument('--targetbam', required=True, help='non-haplotagged bam')
    parser.add_argument('-o', '--outbam', default=None, help='output .bam name (default = (targetbamprefix).copytags.bam)')


    args = parser.parse_args()
    main(args)