## methylartist: transposon methylation nanopore tools

Tools for parsing and plotting methylation data from nanopore data

## Installation

`git clone https://github.com/adamewing/methylartist.git` or download a .zip file from GitHub.

install methylartist and dependencies via:

`python setup.py install`

## Input data

### Alignments (.bam)

Alignments stored in .bam format should be sorted and indexed and should use the same read names as the associated methylation data.

### Methylation Calls

MethylArtist has functions for loading methylation data can from nanopolish, megalodon, or basecalled guppy fast5s, using the appropriate function below.

Once coverted, the sqlite .db file can be input to MethylArtist functions (e.g. `segplot`, `locus`, etc).

## Commands:

### db-nanopolish

Load nanopolish methylation into sqlite db.

Example:

Loading results from `nanopolish call-methylation` to a database:

```
methylartist db-nanopolish -m MCF7_ATCC_REP1.nanopolish.tsv.gz -d MCF7_ATCC.nanopolish.db
```

Appending additional results to the above database:

```
methylartist db-nanopolish -m MCF7_ATCC_REP2.nanopolish.tsv.gz -d MCF7_ATCC.nanopolish.db -a
```

Loading results with the current recommended cutoffs for nanopolish (abs(llr) > 2.0, scale grouped CpGs):

```
methylartist db-nanopolish -m MCF7_ATCC_REP1.nanopolish.tsv.gz,MCF7_ATCC_REP2.nanopolish.tsv.gz -d MCF7_ATCC.nanopolish.db -t 2.0 -s
```

Inputs can be uncompressed or .gzipped.


### db-megalodon

Load megalodon methylation into sqlite db.

The input file is the output of `megalodon_extras per_read_text modified_bases /path/to/megalodon_output`, which needs to be run prior to this script.

The default filename (`/path/to/megalodon_output/per_read_modified_base_calls.txt`) is the same for all megalodon runs, so the `--db` option is recommended to make the output database more identifiable for downstream analysis.

Example:

```
methylartist db-megalodon -m MCF7_ATCC_REP1/per_read_modified_base_calls.txt --db MCF7_ATCC_REP1.db
```

Appending additional results to the above database:

```
methylartist db-megalodon -m MCF7_ATCC_REP2/per_read_modified_base_calls.txt --db MCF7_ATCC_REP2.db
```

Input files can be uncompressed or .gzipped.


### db-guppy

Parses methylation call data from fast5 files where basecalling information has been added by running guppy. Calls are anchored to the reference genome via a user-supplied .bam file aligning the reads derived from the same .fast5 files to a referene genome.

Example:

```
methylartist db-guppy -s MCF7_ATCC -f MCF7_ATCC_REP1/fast5 -p 32 -m [C]G -n 5mC -b MCF7_ATCC.bam -r Homo_sapiens_assembly38.fasta
```

Appending additional results:

```
methylartist db-guppy -s MCF7_ATCC -f MCF7_ATCC_REP2/fast5 -p 32 -m [C]G -n 5mC -b MCF7_ATCC.bam -r Homo_sapiens_assembly38.fasta -a
```

Notes:

The `-f/--fast5` option should point at a folder that contains .fast5 files. Proccessing the fast5 files can be done in parallel with the `-p`/`--procs` option.

Motifs (`-m`/`--motif`) are expressed by putting square brackets around the methylated base e.g. `[C]G` for 5mC and `G[A]TC` for dam-methylated 6mA.

The `-n`/`--modname` option expects the 'long name' stored in the guppy metadata e.g. `5mC` or `6mA`.

Read caching yields a large .db file ending with .bamcache.db. This can be deleted if you to not intend to append additional results to the methylation database (.guppy.db).


### segmeth

Outputs aggregate methylation / demethylation call counts over intervals. Required before generating strip / violin plots with `segplot`

Requires whitespace-delimited list of segments in a BED-like format: chromosome, start, end, label, strand.

Sample input file `-d/--data` has the following whitespace-delimited fields (one sample per line): BAM file, Methylation DB (generated with e.g. `db-nanopolish`)

Highly recommend parallelising with `-p/--proc` option if possible.

Can be used to generate genome-wide methylation stats aggregated over windows via `bedtools makewindows`.

Example:

Aggregate methylation calls over full length L1 elements (RepeatMasker, hg38), exclude reads which map entirely within L1s:

```
methylartist segmeth -d MCF7_data.txt -i L1_FL.bed -p 32 --excl_ambig
```

Contents of `MCF7_data.txt`:

```
MCF7_ATCC.haplotag.bam MCF7_ATCC.nanopolish.db
MCF7_ECACC.haplotag.bam MCF7_ECACC.nanopolish.db
```

Contents of L1_FL.bed (first 10 lines):

```
chr1    440936  447357  L1PA7   +
chr1    675912  682333  L1PA7   +
chr1    7411758 7417922 L1PA4   -
chr1    14355445        14361840        L1PA7   -
chr1    14382319        14388341        L1PA3   -
chr1    25506586        25512707        L1PA2   +
chr1    30427133        30433546        L1MA4   +
chr1    30453586        30459718        L1PA5   -
chr1    30493390        30499556        L1PA4   +
chr1    33671998        33678118        L1PA4   +
```

Output in `L1_FL.MCF7_data.excl_ambig.segmeth.tsv` (first 10 lines):

```
seg_id  seg_chrom       seg_start       seg_end seg_name        seg_strand      MCF7_ATCC.haplotag_meth_calls   MCF7_ATCC.haplotag_unmeth_calls MCF7_ATCC.haplotag_no_calls     MCF7_ATCC.haplotag_methfrac     MCF7_ECACC.haplotag_meth_calls     MCF7_ECACC.haplotag_unmeth_calls MCF7_ECACC.haplotag_no_calls     MCF7_ECACC.haplotag_methfrac
chr1:440936-447357      chr1    440936  447357  L1PA7   +       0       0       0       NaN     0       0       0       NaN
chr1:675912-682333      chr1    675912  682333  L1PA7   +       129     66      141     0.6615384615384615      222     55      133     0.8014440433212996
chr1:7411758-7417922    chr1    7411758 7417922 L1PA4   -       491     410     867     0.5449500554938956      349     79      364     0.8154205607476636
chr1:14355445-14361840  chr1    14355445        14361840        L1PA7   -       222     302     455     0.42366412213740456     293     168     406     0.6355748373101953
chr1:14382319-14388341  chr1    14382319        14388341        L1PA3   -       484     364     557     0.5707547169811321      469     149     563     0.7588996763754046
chr1:25506586-25512707  chr1    25506586        25512707        L1PA2   +       727     495     554     0.5949263502454992      239     276     205     0.4640776699029126
chr1:30427133-30433546  chr1    30427133        30433546        L1MA4   +       35      48      78      0.42168674698795183     20      27      48      0.425531914893617
chr1:30453586-30459718  chr1    30453586        30459718        L1PA5   -       309     297     516     0.5099009900990099      189     42      200     0.8181818181818182
chr1:30493390-30499556  chr1    30493390        30499556        L1PA4   +       460     381     673     0.5469678953626635      328     86      505     0.7922705314009661
```


### segplot

Generates strip plots or violin plots (`-v/--violin`) from `segmeth` output.

Examples:

Strip plot of L1HS, L1PA2, L1PA4, and L1PA6 elements across two samples:

```
methylartist segplot -s L1_FL.MCF7_data.excl_ambig.segmeth.tsv -m MCF7_ATCC.haplotag,MCF7_ECACC.haplotag -c L1HS,L1PA2,L1PA4,L1PA6 --width 10
```

![strip plot](https://github.com/adamewing/methylartist/blob/main/docs/segplot_example_strip.png?raw=true)

As above, but use violin plots:

```
methylartist segplot -s L1_FL.MCF7_data.excl_ambig.segmeth.tsv -m MCF7_ATCC.haplotag,MCF7_ECACC.haplotag -c L1HS,L1PA2,L1PA4,L1PA6 --width 10 --violin
```

![violin plot](https://github.com/adamewing/methylartist/blob/main/docs/segplot_example_violin.png?raw=true)

Note that default output is in .png format. For .svg vector output suitable for editing in inkscape or illustrator add the `--svg` option. Note that for strip plots, this is often inadvisable due to the large number of points.


### locus

Generates smoothed methylation profiles across specific loci with many configurable parameters for one or more samples.

Sample input file `-d/--data` has the following whitespace-delimited fields (one sample per line): BAM file, Methylation DB (generated with e.g. `db-nanopolish`)

Example:

Plot of the GPER1 locus in hg38, highlighting the GeneHancer promoter/enhancer annotation (`GH07J001085`).

```
methylartist locus -d MCF7_data_megalodon.txt -i chr7:1072064-1101499 -g Homo_sapiens.GRCh38.97.chr.sorted.gtf.gz --genes GPER1 -l 1085667-1089471 --cpgspace 30 -p 1,6,1,3,4
```

![locus plot](https://github.com/adamewing/methylartist/blob/main/docs/GPER1.MCF7_data_megalodon.chr7_1072064-1101499.m.locus.meth.png?raw=true)

From top to bottom, the plot shows the genome coordinates, gene models (optional if `--gtf` is supplied), read mappings with modified bases as closed (modified) or open (unmodified) circles, translation from genome to CpG-only coordinate space, raw log-likelihood ratios, and smoothed methylated fraction plot.


### haplocus

Generates smoothed haplotype-aware methylation profiles across specific loci for one sample. Requires haplotypes to be tagged with whatshap.

Examples:

Plot of the PEG3 imprinted region on chromosome 19, hg38.

```
methylartist haplocus -b MCF7_ECACC.haplotag.bam -d MCF7_ECACC.megalodon.db -g Homo_sapiens.GRCh38.97.chr.sorted.gtf.gz --cpgspace 25 -i chr19:56810076-56870725 -l 56835076-56841076 --phasepalette viridis
```

![haplocus plot](https://github.com/adamewing/methylartist/blob/main/docs/MCF7_ECACC.haplotag.chr19_56810076-56870725.m.phased.meth.png?raw=true)


### region

More tractable than "locus" for larger regions using binned methylation.

Example:

Plot of the PGR (Progesterone Receptor) region on chr11 in MCF7 cells. Highlighed region corresponds to the GeneHancer annotation `GH11J101126`. Note use of `-n CG` (requires `-r/--ref` to be set) to normalised bins for content of CpG dinucleotides. Options `--samplepalette` and `--highlightpalette` are used to set colours and `-c 2` has the effect of drawing connectors between genome space and modified base space every 2 bins.

```
methylartist region -i chr11:100956822-101228191 -d MCF7_data_megalodon.txt -g Homo_sapiens.GRCh38.97.chr.sorted.gtf.gz -p 32 -n CG -r /home/data/ref/hg38/Homo_sapiens_assembly38.fasta --samplepalette magma -l 101126888-101129371 --highlightpalette viridis -c 2
```

![region plot](https://github.com/adamewing/methylartist/blob/main/docs/MCF7_data_megalodon.chr11_100956822-101228191.m.region.meth.png?raw=true)


### composite

Generates "composite" methylation plots where multiple per-element profiles are aligned to and plotted against a reference sequence.

Example:

First aggregate methylation counts over intervals corresponding to the TE family of interest:
```
methylartist segmeth -d MCF7_data_megalodon.txt -i L1HS.bed -p 32 --excl_ambig
```

The contents of `L1HS.bed` are formatted as follows (first 10 lines), derived from hg38 repeatmasker .out files:

```
chr1    34566056        34572105        L1HS    -
chr1    67078892        67084915        L1HS    -
chr1    71513699        71519742        L1HS    +
chr1    80939204        80945257        L1HS    -
chr1    84052390        84058406        L1HS    -
chr1    85748520        85754548        L1HS    +
chr1    85927068        85933100        L1HS    +
chr1    86679081        86685111        L1HS    -
chr1    104770248       104776278       L1HS    -
chr1    104843834       104849864       L1HS    -
```

(adjust number of processes used for computing individual profiles via `-p/--procs` to be appropriate for your system)

The consensus sequence for establishing a common coordinate system is in `L1.3.fa` (fasta formatted with one sequence), available here: [L19088.1](https://www.ncbi.nlm.nih.gov/nuccore/L19088.1)

The contents of `L1.3.highlights.bed` used for optional highlighting of regions (ORF1 and ORF2 in this case):
```
910 1926 ORF1 #bcbddc
1990 5817 ORF2 #756bb1
```

Output file from `methylartist segmeth`, `L1HS.MCF7_data_megalodon.excl_ambig.segmeth.tsv` (first 10 lines):

```
seg_id  seg_chrom       seg_start       seg_end seg_name        seg_strand      MCF7_ATCC.haplotag_m_meth_calls MCF7_ATCC.haplotag_m_unmeth_calls       MCF7_ATCC.haplotag_m_no_calls   MCF7_ATCC.haplotag_m_methfrac    MCF7_ECACC.haplotag_m_meth_calls        MCF7_ECACC.haplotag_m_unmeth_calls      MCF7_ECACC.haplotag_m_no_calls  MCF7_ECACC.haplotag_m_methfrac
chr1:34566056-34572105  chr1    34566056        34572105        L1HS    -       1121    1182    399     0.4867564046895354      891     613     276     0.5924202127659575
chr1:67078892-67084915  chr1    67078892        67084915        L1HS    -       1244    726     303     0.6314720812182741      1229    686     287     0.64177545691906
chr1:71513699-71519742  chr1    71513699        71519742        L1HS    +       881     687     340     0.5618622448979592      760     509     248     0.598896769109535
chr1:80939204-80945257  chr1    80939204        80945257        L1HS    -       87      19      13      0.8207547169811321      669     258     198     0.7216828478964401
chr1:84052390-84058406  chr1    84052390        84058406        L1HS    -       223     103     46      0.6840490797546013      0       0       0       NaN
chr1:85748520-85754548  chr1    85748520        85754548        L1HS    +       902     503     265     0.6419928825622776      861     482     212     0.6411020104244229
chr1:85927068-85933100  chr1    85927068        85933100        L1HS    +       1407    630     265     0.6907216494845361      443     288     118     0.6060191518467852
chr1:86679081-86685111  chr1    86679081        86685111        L1HS    -       717     468     223     0.6050632911392405      758     386     185     0.6625874125874126
chr1:104770248-104776278        chr1    104770248       104776278       L1HS    -       546     664     248     0.4512396694214876      828     829     256     0.4996982498491249
```

Command to `methylartist composite` (note the parsing of individual profiles can be parallelised via `-p`):

```
methylartist composite -b MCF7_ATCC.haplotag.bam -m MCF7_ATCC.megalodon.db --sample MCF7_ATCC.haplotag_m -s L1HS.MCF7_data_megalodon.excl_ambig.segmeth.tsv -f L1HS -r Homo_sapiens_assembly38.fasta -t L1.3.fa -p 32 --blocks L1.3.highlights.bed
```

![composite plot](https://github.com/adamewing/methylartist/blob/main/docs/MCF7_ATCC.haplotag.L1HS.composite.png?raw=true)


### wgmeth

Outputs genome-wide statistics on CpGs covered by at least one call. By default, `wgmeth` yields output in bedMethyl format (0-based coordinates):
```
chr1    10467   10468   MCF7_ECACC.nanopolish    3       .       10467   10468   0,0,0   3       1.0
chr1    10469   10470   MCF7_ECACC.nanopolish    3       .       10469   10470   0,0,0   3       1.0
chr1    10482   10483   MCF7_ECACC.nanopolish    7       .       10482   10483   0,0,0   7       0.7142857142857143
chr1    10487   10488   MCF7_ECACC.nanopolish    7       .       10487   10488   0,0,0   7       0.7142857142857143
chr1    10491   10492   MCF7_ECACC.nanopolish    7       .       10491   10492   0,0,0   7       0.7142857142857143
chr1    10495   10496   MCF7_ECACC.nanopolish    7       .       10495   10496   0,0,0   7       0.7142857142857143
chr1    10523   10524   MCF7_ECACC.nanopolish    8       .       10523   10524   0,0,0   8       0.875
chr1    10540   10541   MCF7_ECACC.nanopolish    11      .       10540   10541   0,0,0   11      0.9090909090909091
chr1    10561   10562   MCF7_ECACC.nanopolish    10      .       10561   10562   0,0,0   10      1.0
chr1    10569   10570   MCF7_ECACC.nanopolish    10      .       10569   10570   0,0,0   10      1.0
```

The `--dss` option yields genome-wide files suitable for input to DSS for statistical assessment of differential methylation (1-based coordinates):

```
chr     pos     N       X
chr1    10468   3       3
chr1    10470   3       3
chr1    10483   6       5
chr1    10488   6       5
chr1    10492   6       5
chr1    10496   6       5
chr1    10524   6       6
chr1    10541   8       7
chr1    10562   8       8
```
