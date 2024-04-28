#!/usr/bin/env python

import os
import sys
import argparse
import logging
import multiprocessing as mp
import sqlite3

from itertools import product, chain
from operator import itemgetter
from collections import defaultdict as dd

import numpy as np
import scipy.spatial.distance as distance
import networkx as nx
import pysam

FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def overlap(a, b):
    assert len(a) == len(b) == 2
    return min(max(a), max(b)) - max(min(a), min(b))


def smooth(seg, hap, x, window_len=8, window='hanning'):
    # used for locus plots
    ''' modified from scipy cookbook: https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html '''

    assert window_len % 2 == 0, '--smoothwindowsize must be an even number'
    assert x.ndim == 1

    if x.size <= window_len:
        logger.warning('%s: hap %s number of data points not greater than window length' % (seg, hap))
        return x

    if window_len < 3:
        return x

    assert window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')

    return y[(int(window_len/2)-1):-(int(window_len/2))]


def get_modnames(meth_db):
    bam = pysam.AlignmentFile(meth_db)
    logger.info(f'searching for mod names in {meth_db}, if this takes a long time please ensure MM/ML tags are present')

    mod_names = []

    for read in bam.fetch():
        if read.modified_bases is not None:
            for k in read.modified_bases.keys():
                mod_names.append(k[2])

        if len(mod_names) > 0:
            break

    return mod_names


def map_hap_ends(hap_reads, mod_sites, chrom, start, end, min_ext_count=200):
    seg = '%s:%d-%d' % (chrom, start, end)
    assert len(hap_reads) > 2, '%s: hap_reads <= 2' % seg

    L = []
    R = []

    read_bounds = {}
    hap_bounds = {}
    min_ext = {}
    skip = {}

    for hap in hap_reads:
        hap_b = []
        read_b = []

        for readname in hap_reads[hap]:
            read_b += hap_reads[hap][readname]
            hap_b += map(int, list(mod_sites[readname].keys()))

        if len(hap_b) == 0:
            continue

        read_bounds[hap] = [min(read_b), max(read_b)]
        hap_bounds[hap] = [min(hap_b), max(hap_b)]

        hap_b.sort()
        
        ext_count = min_ext_count

        if len(hap_b) < min_ext_count:
            logger.warning('%s: short haplotype: %s (%d called sites)' % (seg, hap, len(hap_b)))
            ext_count = len(hap_b)-1

        if len(hap_b) < 50: # TODO param
            skip[hap] = len(hap_b)

        min_ext[hap] = [hap_b[ext_count]-hap_bounds[hap][0], hap_bounds[hap][1]-hap_b[-ext_count]]

    hap_order = sorted([[hap, hap_bounds[hap][0]] for hap in hap_bounds], key=itemgetter(1))
    #print(hap_order)

    next_hap_PS = {}
    for i, hap in enumerate(hap_order[1:],1):
        prev_hap_PS = hap_order[i-1][0].split(':')[1]
        this_hap_PS = hap_order[i][0].split(':')[1]

        if prev_hap_PS != this_hap_PS:
            next_hap_PS[prev_hap_PS] = this_hap_PS # would overwrite previous assignment (intended)

    #print(next_hap_PS)
    
    #print(min_ext)

    for hap in hap_bounds:
        if hap in skip:
            logger.warning('%s: skip very short haplotype: %s (%d calls)' % (seg, hap, skip[hap]))
            continue

        if read_bounds[hap][0] > start and read_bounds[hap][0] < end:
            L.append(hap)
        
        if read_bounds[hap][1] > start and read_bounds[hap][1] < end:
            R.append(hap)

    forward_join = {}

    for end_l, end_r in product(L,R):
        end_l_PS = end_l.split(':')[1]
        end_r_PS = end_r.split(':')[1]

        if end_r_PS in next_hap_PS:
            if next_hap_PS[end_r_PS] == end_l_PS:
                ovl = overlap(hap_bounds[end_l], hap_bounds[end_r])
                extent = max(ovl, min_ext[end_l][0], min_ext[end_r][1])

                meth_r_bounds = [hap_bounds[end_r][1]-extent,hap_bounds[end_r][1]]
                meth_l_bounds = [hap_bounds[end_l][0],hap_bounds[end_l][0]+extent]

                logger.debug('%s (R:%s, L:%s): extent = %d  overlap = %d' % (seg, end_r, end_l, extent, ovl))
                logger.debug('%s (R:%s, L:%s): meth_r_bounds: %d-%d' % (seg, end_r, end_l, meth_r_bounds[0], meth_r_bounds[1]))
                logger.debug('%s (R:%s, L:%s): meth_l_bounds: %d-%d' % (seg, end_r, end_l, meth_l_bounds[0], meth_l_bounds[1]))

                profile_r = meth_profile(seg, end_r, hap_reads, mod_sites, meth_r_bounds[0], meth_r_bounds[1])
                profile_l = meth_profile(seg, end_l, hap_reads, mod_sites, meth_l_bounds[0], meth_l_bounds[1])

                if len(profile_l) > len(profile_r):
                    profile_l = profile_l[-len(profile_r):]

                elif len(profile_r) > len(profile_l):
                    profile_r = profile_r[:len(profile_l)]

                d_euc = distance.euclidean(profile_l, profile_r)

                if end_r not in forward_join:
                    forward_join[end_r] = (end_l, d_euc)
                
                else:
                    if forward_join[end_r][1] > d_euc:
                        forward_join[end_r] = (end_l, d_euc)

    # remove repeat left ends
    reverse_join = {}

    for hap in forward_join:
        end_l, d_euc = forward_join[hap]

        if end_l not in reverse_join:
            reverse_join[end_l] = (hap, d_euc)
        
        else:
            if reverse_join[end_l][1] > d_euc:
                reverse_join[end_l] = (hap, d_euc)

    forward_join = {}

    for hap in reverse_join:
        end_r, d_euc = reverse_join[hap]

        forward_join[end_r] = (hap, d_euc)

    # if PS:1 --> PS_next:2 then PS:2 --> PS_next:1

    new_forward_join = {}

    for hap in forward_join:
        new_forward_join[hap] = forward_join[hap]
        HP, PS = hap.split(':')
        HP_next, PS_next = forward_join[hap][0].split(':')

        assert HP in ('1','2')
        assert HP_next in ('1','2')

        other_HP = '1'
        other_HP_next = '1'

        if HP == '1':
            other_HP = '2'
        
        if HP_next == '1':
            other_HP_next = '2'

        other_hap = '%s:%s' % (other_HP, PS)
        other_hap_next = '%s:%s' % (other_HP_next, PS_next)

        if other_hap not in forward_join and other_hap in min_ext and other_hap_next in min_ext:
            end_r = other_hap
            end_l = other_hap_next

            ovl = overlap(hap_bounds[end_l], hap_bounds[end_r])
            extent = max(ovl, min_ext[end_l][0], min_ext[end_r][1])

            meth_r_bounds = [hap_bounds[end_r][1]-extent,hap_bounds[end_r][1]]
            meth_l_bounds = [hap_bounds[end_l][0],hap_bounds[end_l][0]+extent]

            profile_l = meth_profile(seg, end_l, hap_reads, mod_sites, meth_l_bounds[0], meth_l_bounds[1])
            profile_r = meth_profile(seg, end_r, hap_reads, mod_sites, meth_r_bounds[0], meth_r_bounds[1])

            if len(profile_l) > len(profile_r):
                profile_l = profile_l[-len(profile_r):]

            elif len(profile_r) > len(profile_l):
                profile_r = profile_r[:len(profile_l)]

            d_euc = distance.euclidean(profile_l, profile_r)

            new_forward_join[end_r] = (end_l, d_euc)

    logger.info('%s: missing neighbor search found %d pairs' % (seg, len(new_forward_join)-len(forward_join)))

    forward_join = new_forward_join

    hapmap = []

    for hap in forward_join:
        link = (hap, forward_join[hap][0], forward_join[hap][1])
        hapmap.append(link)

    hapgraph = nx.Graph()
    hapgraph.add_weighted_edges_from(hapmap)
    deg = hapgraph.degree()

    # prune nodes with degree > 2
    for hap in hapgraph.nodes():
        if deg[hap] > 2:
            logger.warning('%s: node %s has degree %d' % (seg, hap, deg[hap]))
            edgelist = []
            for n_hap in hapgraph.neighbors(hap):
                w = float(hapgraph[hap][n_hap]['weight'])
                edgelist.append((hap, n_hap, w))
            
            for e in sorted(edgelist, key=itemgetter(2))[2:]:
                hapgraph.remove_edge(e[0], e[1])
                logger.warning('%s: removed edge %s --> %s (%.2f)' % (seg, e[0], e[1], e[2]))

    try:
        for cyc in nx.algorithms.cycles.find_cycle(hapgraph):
            print('%s cycle in haplotype graph! debug:', cyc)
    except nx.exception.NetworkXNoCycle:
        logger.info('%s: no cycles in haplotype graph' % seg)
        
    ccs = list(nx.connected_components(hapgraph))

    logger.info('%s: initial haplotype graph has %d nodes' % (seg, hapgraph.number_of_nodes()))
    logger.info('%s: initial haplotype graph has %d connected components' % (seg, len(ccs)))


    # close gaps based on proximity (generally non-overlapping)

    iter = 0
    
    while len(ccs) > 1: # run on len(ccs) == 2 for stats

        incompat_haps = []

        for i, cc in enumerate(ccs):
            for hap in cc:
                HP, PS = hap.split(':')

                other_HP = '1'
                if HP == '1':
                    other_HP = '2'
                
                other_hap = '%s:%s' % (other_HP, PS)

                if other_hap in cc:
                    incompat_haps.append((hap, other_hap, i))

        logger.info('%s: %d self-incompatible haplotypes in initial components' % (seg, len(incompat_haps)))

        compat_ccs = {}

        r = list(range(len(ccs)))
        checked_r = []

        for i,j in product(r,r):
            if i == j:
                continue

            pair = ','.join(map(str, sorted((i,j))))

            if pair in checked_r:
                continue

            checked_r.append(pair)
            compat = check_compat_ccs(ccs[i], ccs[j])
            compat_ccs[pair] = compat

        end_coords = []

        for cc in ccs:
            c_end_coords = []

            for c in cc:
                assert deg[c] <= 2

                if(deg[c]) <= 1:
                    c_end_coords += hap_bounds[c]
            
            end_coords.append([min(c_end_coords), max(c_end_coords)])

        assert len(end_coords) == len(ccs)

        compat_overlaps = {}
        contained_ccs = dd(dict)

        for i, j in product(r,r):
            if i == j:
                continue

            pair = ','.join(map(str, sorted((i,j))))

            if pair in compat_ccs and compat_ccs[pair] == 0:
                compat_overlaps[pair] = overlap(end_coords[i], end_coords[j])

        # find the ends of each haplotype chain

        min_2_coords = []
        max_2_coords = []

        min_2_index = []
        max_2_index = []

        for i, ec in enumerate(end_coords):
            if len(min_2_coords) < 2:
                min_2_coords.append(ec[0])
                max_2_coords.append(ec[1])
                min_2_index.append(i)
                max_2_index.append(i)
                continue
            
            if ec[0] < min_2_coords[0]:
                min_2_coords[0] = ec[0]
                min_2_index[0] = i
            elif ec[0] < min_2_coords[1]:
                min_2_coords[1] = ec[0]
                min_2_index[1] = i

            if ec[1] > max_2_coords[0]:
                max_2_coords[0] = ec[1]
                max_2_index[0] = i
            elif ec[1] > max_2_coords[1]:
                max_2_coords[1] = ec[1]
                max_2_index[1] = i

        min_2_incompat = compat_ccs[','.join(map(str, sorted(min_2_index)))]
        max_2_incompat = compat_ccs[','.join(map(str, sorted(max_2_index)))]

        print(compat_ccs)
        print(end_coords)
        print(compat_overlaps)
        print(min_2_coords, min_2_index, min_2_incompat)
        print(max_2_coords, max_2_index, max_2_incompat)

        logger.info('%s: minimum joined haplotype end coords: %d, %d' % (seg, min_2_coords[0], min_2_coords[1]))
        logger.info('%s: maximum joined haplotype end coords: %d, %d' % (seg, max_2_coords[0], max_2_coords[1]))


        if len(ccs) > 2:
            best_join = None

            for pair_str, dist in compat_overlaps.items():
                dist = abs(dist)
                pair = list(map(int, pair_str.split(',')))

                # if pair[1] in min_2_index:
                #     continue

                # if pair[0] in max_2_index:
                #     continue

                if best_join is None:
                    best_join = pair

                else:
                    best_pair_str = ','.join(map(str, sorted(best_join)))
                    if dist < abs(compat_overlaps[best_pair_str]):
                        best_join = pair
            
            end_r = None
            end_l = None

            if best_join is None:
                logger.warning('%s: best_join is None' % seg)
                print('compat_overlaps:', compat_overlaps)
                print('min_2_index:', min_2_index)
                print('min_2_coords', min_2_coords)
                print('max_2_index:', max_2_index)
                print('max_2_coords', max_2_coords)
                print('end_coords', end_coords)

                sys.exit(1)

            for c in ccs[best_join[0]]:
                if hap_bounds[c][1] == end_coords[best_join[0]][1]:
                    end_r = c

            for c in ccs[best_join[1]]:
                if hap_bounds[c][0] == end_coords[best_join[1]][0]:
                    end_l = c

            pos_r = hap_bounds[end_r][1]
            pos_l = hap_bounds[end_l][0]

            logger.info('%s: best join: %s (%d) --> %s (%d)' % (seg, end_r, pos_r, end_l, pos_l))

            link = (end_r, end_l, 99.0)
            hapgraph.add_weighted_edges_from([link])

            ccs = list(nx.connected_components(hapgraph))
        
            logger.info('%s: haplotype graph has %d connected components after iteration %d' % (seg, len(ccs), iter))

            iter += 1

        else:
            break


    # if hap 1:xxx in cc 0 but hap2:xxx not in cc 1 put hap2:xxx in cc 1 and vice versa

    assert len(ccs) == 2

    add_cc = {}

    for i, cc in enumerate(ccs):
        for hap in cc:
            HP, PS = hap.split(':')

            other_HP = '1'
            if HP == '1':
                other_HP = '2'

            j = 0
            if i == 0:
                j = 1
            
            other_cc = ccs[j]
            
            other_hap = '%s:%s' % (other_HP, PS)

            assert other_hap not in cc

            if other_hap not in other_cc:
                logger.info('%s: %s in cc %d: adding opposite hap %s to other cc' % (seg, hap, i, other_hap))

                ccs[j] = other_cc | set([other_hap])

    cc_remap = {}

    for i, cc in enumerate(ccs):
        HP = str(i+1)

        for hap in cc:
            cc_remap[hap] = int(HP)

    return chrom, cc_remap


def check_compat_ccs(cc1, cc2):
    ic = 0

    for hap in cc1:
        HP, PS = hap.split(':')

        other_HP = '1'
        if HP == '1':
            other_HP = '2'
        
        other_hap = '%s:%s' % (other_HP, PS)

        if other_hap in cc2:
            ic += 1
    
    return ic
        

def meth_profile(seg, hap, hap_reads, mod_sites, start, end):
    meth_table = dd(dict) # site-->read-->call

    for readname in hap_reads[hap]:
        for site in mod_sites[readname]:
            if int(site) >= start and int(site) < end:
                meth_call = mod_sites[readname][site]
                if meth_call in (-1, 1):
                    meth_table[int(site)][readname] = meth_call
    
    raw_profile = []
    for pos in sorted(list(meth_table.keys())):
        p = list(meth_table[pos].values())
        raw_profile.append(p.count(1)/(p.count(1) + p.count(-1)))

    sm_window = int(len(raw_profile)/100)
    if sm_window % 2 != 0:
        sm_window += 1 # must be even

    if sm_window < 8:
        sm_window = 8

    sm_profile = smooth(seg, hap, np.asarray(raw_profile), window_len=sm_window)

    return sm_profile


def hap_read_dict(bam_fn, mod, chrom, start, end, pad=10000):
    bam = pysam.AlignmentFile(bam_fn)

    start = int(start + pad)
    end = int(end - pad)

    haps = dd(dict)
    mod_sites = dd(dict)

    for read in bam.fetch(chrom, start, end):
        HP = None
        PS = None

        if read.is_supplementary or read.is_secondary or read.is_duplicate:
            continue

        if read.modified_bases is None:
            continue

        for tag in read.get_tags():
            if tag[0] == 'HP':
                HP = tag[1]
            if tag[0] == 'PS':
                PS = tag[1]
            
        if None not in (HP,PS):
            if HP not in (1,2):
                logger.warning('found HP %d, this is designed for diploid haplotypes only' % HP)
            hap = '%d:%d' % (HP,PS)

            if hap in haps:
                if read.query_name in haps[hap]:
                    if haps[hap][read.query_name][0] > read.reference_start:
                        haps[hap][read.query_name][0] = read.reference_start

                    if haps[hap][read.query_name][1] < read.reference_end:
                        haps[hap][read.query_name][1] = read.reference_end

                else:
                    haps[hap][read.query_name] = [read.reference_start, read.reference_end]

            else:
                haps[hap][read.query_name] = [read.reference_start, read.reference_end]

        for modtype in read.modified_bases:
            if mod != modtype[2]:
                continue

            ref_map = dict(read.get_aligned_pairs(matches_only=True))

            for (read_bp, score) in read.modified_bases[modtype]:
                if read_bp not in ref_map:
                    continue

                ref_bp = ref_map[read_bp]

                if ref_bp < start or ref_bp > end:
                    continue

                methstate = 0
                if score/255 > 0.8:
                    methstate = 1

                if score/255 < 0.2:
                    methstate = -1

                mod_sites[read.query_name][ref_bp] = methstate

    return haps, mod_sites


def match_hap_ends(args, chrom, start, end):
    hap_reads, mod_sites = hap_read_dict(args.bam, args.mod, chrom, start, end)

    if len(hap_reads) > 2:
        return map_hap_ends(hap_reads, mod_sites, chrom, start, end)

    return None, None



def rehap_bam(args, hapmap):
    out_fn = '.'.join(args.bam.split('.')[:-1]) + '.methmergehaps.bam'

    logger.info('writing reads to %s' % out_fn)

    in_bam = pysam.AlignmentFile(args.bam)
    out_bam = pysam.AlignmentFile(out_fn, 'wb', template=in_bam)

    change_count = 0
    read_count = 0

    for read in in_bam.fetch(until_eof=True):
        HP = None
        PS = None

        for tag in read.get_tags():
            if tag[0] == 'HP':
                HP = tag[1]
            if tag[0] == 'PS':
                PS = tag[1]
            
        if None not in (HP,PS):
            if HP not in (1,2):
                logger.warning('found HP %d this is designed for diploid haplotypes only' % HP)
            hap = '%d:%d' % (HP,PS)

            if read.reference_name in hapmap:
                if hap in hapmap[read.reference_name]:
                    read.set_tag('HP', hapmap[read.reference_name][hap], value_type='i')
                    if HP != hapmap[read.reference_name][hap]:
                        change_count += 1

        out_bam.write(read)
        read_count += 1

    out_bam.close()
    logger.info('wrote %d reads to %s, changed %d HP tags' % (read_count, out_fn, change_count))


def main(args):
    avail_mods = get_modnames(args.bam)
    
    if args.mod not in avail_mods:
        logger.error('--mod %s not available. Possible options: %s' % (args.mod, ','.join(list(set(avail_mods)))))
        sys.exit(1)
    
    chromlen = {}
    with open(args.fai) as fai:
        for line in fai:
            chrom, length = line.strip().split()[:2]
            length = int(length)
            if length > int(args.minchromlen):
                chromlen[chrom] = length

    logger.info('%d chromosomes with length > %d' % (len(chromlen), int(args.minchromlen)))

    if args.debug_seg:
        chrom, span = args.debug_seg.split(':')
        start, end = map(int, span.split('-'))

        result = match_hap_ends(args, chrom, start, end)

        sys.exit()

    pool = mp.Pool(processes=int(args.procs))

    results = []

    for chrom in chromlen:
        res = pool.apply_async(match_hap_ends, [args, chrom, 0, chromlen[chrom]])
        results.append(res)

    hapmap = {}

    for res in results:
        chrom, cc_remap = res.get()
        
        if chrom is not None:
            hapmap[chrom] = cc_remap

    rehap_bam(args, hapmap)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='methylartist: tools for exploring nanopore modified base data')
    parser.add_argument('-b', '--bam', required=True, help='indexed sorted bam file with MM/ML tags')
    parser.add_argument('-m', '--mod', required=True, help='modificaiton (if not sure, guess and then select from list if error)')
    parser.add_argument('-f', '--fai', required=True, help='.fai file from "samtools faidx" or equivalent')
    parser.add_argument('-p', '--procs', default=4, help='number of processes (chromosomes) to run concurrently (default = 4)')
    parser.add_argument('--debug_seg', default=None)
    parser.add_argument('--minchromlen', default='10000000', help='minimum chromosome size for consideration (default=10000000)')
    

    args = parser.parse_args()
    main(args)
