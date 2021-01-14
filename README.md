## tmnt: Transposon methylation nanopore tools

Tools for parsing and plotting methylation data from Oxford Nanopore Technologies data in a way that is useful for transposable element applications.

## Installation

`git clone https://github.com/adamewing/tmnt.git` or download a .zip file from GitHub.

install tmnt and dependencies via:

`python setup.py install`

## Input data

### Alignments (.bam)

Alignments stored in .bam format should be sorted and indexed and should use the same read names as the associated methylation data.

### Methylation Calls

TMNT hs functions for loading methylation data can from nanopolish, megalodon, or basecalled guppy fast5s, using the appropriate function below.

Once coverted, the sqlite .db file can be input to TMNT functions (e.g. `segplot`, `locus`, etc).

## Commands:

### db-nanopolish

Load nanopolish methylation into sqlite db. Needs to be done before running methylation parsing or plotting commands.

### db-megalodon
(work in progress)

### db-guppy
(work in progress)

### segmeth

Outputs aggregate methylation / demethylation call counts over intervals. Required before generating strip / violin plots with `segplot`

Requires whitespace-delimited list of segments in a BED-like format: chromosome, start, end, label, strand.

Sample input file `-d/--data` has the following whitespace-delimited fields (one sample per line): BAM file, Methylation DB (generated with e.g. `db-nanopolish`)

Highly recommend parallelising with `-p/--proc` option if possible.

### segplot

Generates strip plots or violin plots (`-v/--violin`) from `segmeth` output.

### locus

Generates smoothed methylation profiles across specific loci with many configurable parameters for one or more samples.

Sample input file `-d/--data` has the following whitespace-delimited fields (one sample per line): BAM file, Methylation DB (generated with e.g. `db-nanopolish`)

### haplocus

Generates smoothed haplotype-aware methylation profiles across specific loci for one sample. Requires haplotypes to be tagged with whatshap.

### composite

Generates "composite" methylation plots where multiple per-element profiles are aligned to and plotted against a reference sequence.

### dss

Creates genome-wide files suitable for input to DSS for statistical assessment of differential methylation.
