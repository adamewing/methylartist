#!/usr/bin/env python

import sys
import logging
import argparse
import pandas as pd
import scipy.stats as ss
import scikit_posthocs as sp
import numpy as np

from itertools import combinations

FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def single_arg(args, f):
    data = pd.read_csv(args.tsv, sep='\t', header=0)

    avail_annots = set(data['seg_name'])
    annots = list(avail_annots)

    if args.annotations:
        annots = args.annotations.split(',')
        for ann in annots:
            if ann not in avail_annots:
                logger.error(f'{ann} not found in annotations. Available annotations:')
                logger.error('\n'+'\n'.join(list(avail_annots)))
                sys.exit(1)
    
    cols = data.columns[6:]

    if args.columns:
        cols = args.columns.split(',')
        for col in cols:
            if col not in data.columns:
                logger.error(f'{col} not found in columns. Available columns:')
                logger.error('\n' + '\n'.join(data.columns))
                sys.exit(1)

    for ann in annots:
        subdata = data[data['seg_name']==ann]
        subdata = subdata[cols]
        res = f(subdata)

        print(f'\n{ann}:')

        with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            print(res.to_string())


def mean(args):
    single_arg(args, pd.DataFrame.mean)


def median(args):
    single_arg(args, pd.DataFrame.median)


def kruskal(args):
    data = pd.read_csv(args.tsv, sep='\t', header=0)

    avail_annots = set(data['seg_name'])
    annots = list(avail_annots)

    if args.annotations:
        annots = args.annotations.split(',')
        for ann in annots:
            if ann not in avail_annots:
                logger.error(f'{ann} not found in annotations. Available annotations:')
                logger.error('\n'+'\n'.join(list(avail_annots)))
                sys.exit(1)

    cols = args.columns.split(',')

    for col in cols:
        if col not in data.columns:
            logger.error(f'{col} not found in columns. Available columns:')
            logger.error('\n'+'\n'.join(data.columns))
            sys.exit(1)
    
    if len(cols) < 2:
        sys.exit('kruskal requires at least two columns')

    if not (args.posthoc):
        print('Annotation\tH_Stat\tP_Value')

    for ann in annots:
        subdata = data[data['seg_name']==ann]
        subdata = subdata[cols]
        h_stat, p_value = ss.kruskal(*[subdata[col].dropna() for col in subdata.columns])

        if args.posthoc:
            print('Annotation\tH_Stat\tP_Value')


        print(f'{ann}\t{h_stat}\t{p_value}')

        if args.posthoc:
            print(f'\nPost-hoc Dunn\'s tests:')
            sd_long = pd.melt(subdata, var_name='group', value_name='values')
            posthoc = sp.posthoc_dunn(sd_long, val_col='values', group_col='group', p_adjust='holm')
            print(posthoc)
            print('\n')


def effsize(args):
    data = pd.read_csv(args.tsv, sep='\t', header=0)

    avail_annots = set(data['seg_name'])
    annots = list(avail_annots)

    if args.annotations:
        annots = args.annotations.split(',')
        for ann in annots:
            if ann not in avail_annots:
                logger.error(f'{ann} not found in annotations. Available annotations:')
                logger.error('\n'+'\n'.join(list(avail_annots)))
                sys.exit(1)

    cols = args.columns.split(',')

    for col in cols:
        if col not in data.columns:
            logger.error(f'{col} not found in columns. Available columns:')
            logger.error('\n'+'\n'.join(data.columns))
            sys.exit(1)
    
    if len(cols) < 2:
        sys.exit('must specify at least two columns')
    
    for ann in annots:
        subdata = data[data['seg_name']==ann]
        subdata = subdata[cols]

        for c1, c2 in combinations(cols, 2):
            d = cohen_d(subdata[c1], subdata[c2])
            d = '%.3f' % d
            print(f'{ann}\t{c1}\t{c2}\t{d}')
    

def cohen_d(m1, m2):
    return (np.mean(m1) - np.mean(m2))/np.var(pd.concat([m1, m2]))


def main():
    logger.info('segstats cmd: %s' % ' '.join(sys.argv))
    args = parse_args()
    args.func(args)


def parse_args():
    parser = argparse.ArgumentParser(description='simple stats on segplot output')

    parser.add_argument('tsv', help='segmeth .tsv')
    parser.add_argument('-a', '--annotations', default=None, help='limit analysis to one or more (comma-delimited) annotations (seg_name column)')

    subparsers = parser.add_subparsers(title="command", dest="command")
    subparsers.required = True

    parser_mean = subparsers.add_parser('mean')
    parser_median = subparsers.add_parser('median')
    parser_kruskal = subparsers.add_parser('kruskal')
    parser_effsize = subparsers.add_parser('effsize')

    parser_mean.set_defaults(func=mean)
    parser_median.set_defaults(func=median)
    parser_kruskal.set_defaults(func=kruskal)
    parser_effsize.set_defaults(func=effsize)

    parser_mean.add_argument('-c', '--columns', default=None, help='which columns to output (comma-delimited, default = all)')
    parser_median.add_argument('-c', '--columns', default=None, help='which columns to output (comma-delimited, default = all)')
    parser_kruskal.add_argument('-c', '--columns', required=True, help='comma-delimited columns')
    parser_kruskal.add_argument('--posthoc', action='store_true', default=False, help='perform post-hoc Dunn\'s tests')
    parser_effsize.add_argument('-c', '--columns', required=True, help='columns to compare (pairwise comparisons)')

    return parser.parse_args()

if __name__ == '__main__':
    main()