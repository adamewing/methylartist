.. _commands:

Commands
=========

The ``methylartist`` command-line suite provides a collection of tools for exploring nanopore-derived modified base data. It is built with Python and leverages packages such as:

  - **pysam** – for parsing sequence alignment map (BAM) and methylation database files  
  - **matplotlib** and **seaborn** – for generating detailed plots of methylation distribution  

``methylartist`` works with standard biological file formats including **BAM**, gzipped **VCF**, and **GTF**.  
All input files must be properly indexed and index files named according to the original filename.  

.. note::
   For GTF files, the chromosome column should use the format ``chr10`` rather than ``1``.

To view all available subcommands, run:

.. code-block:: console

   methylartist -h

You should see output similar to the following:

.. code-block:: console

   options:
     -h, --help            show this help message and exit
     -v, --version         show program's version number and exit

   tool:
     {db-nanopolish,db-megalodon,db-custom,db-guppy,db-sub,
      segmeth,segplot,locus,region,composite,wgmeth,
      adjustcutoffs,scoredist}

``methylartist`` currently includes **13 subcommands**. Each subcommand is designed for specific tasks, ranging from database creation to visualisation and analysis of methylation data. 

.. note::
   The five subcommands prefixed with ``db-`` create methylation databases from different input formats. These databases are required for visualisation when BAM files lack ``MM`` and ``ML`` tags.

The sections below describe each subcommand in detail, including its purpose, parameters, and example usage. These parameters can be viewed using ``subcommand -h``.  

db-nanopolish
---------------

The ``db-nanopolish`` subcommand loads **nanopolish methylation data** into an SQLite database.

**Options:**

.. option:: -m, --methdata 
  
   Whole-genome Nanopolish methylation output file. Multiple files can be provided as a comma-delimited list.

.. option:: -d, --db  
   
   Preferred name of the output SQLite database (the default is auto-inferred from input).

.. option:: -t, --thresh

   Log-likelihood ratio (LLR) threshold for methylation calls. The default is 2.5.  
   If using --scalegroup, the suggested setting is 2.0.

.. option:: -a, --append

   Append results to an existing database (.db) file.

.. option:: -s, --scalegroup

   Scale threshold by the number of CpGs in a group.

.. option:: -n, --modname

   Modification name or tag (useful when combining multiple modification types). Default is CpG.

.. option:: --motif

   Methylation motif to use. Default is CG.


**Example usage:**

.. code-block:: console

   methylartist db-nanopolish -m sample_methylation.tsv -d sample_methyl.db

This command creates a methylation database named ``sample_methyl.db`` from the provided Nanopolish output.



db-megalodon
-------------

The ``db-megalodon`` subcommand loads **Megalodon methylation output** into an SQLite database.  
This command processes the output from the following Megalodon command:

.. code-block:: console

   megalodon_extras per_read_text modified_bases /path/to/megalodon_output

**Options:**

.. option:: -m, --methdata

   Megalodon per_read_text methylation output file.

.. option:: -d, --db

   Name of the output SQLite database (default is auto-inferred from input).

.. option:: -p, --minprob

   Probability threshold for calling a base as modified or unmodified. Default is 0.8.

.. option:: -a, --append

   Append results to an existing database.

.. option:: --motifsize

   Motif size for the modified base. Default is 2 (appropriate for the common CG context).  
   Set to 1 for single-base modifications such as 6mA.



**Example usage:**

.. code-block:: console

   methylartist db-megalodon -m sample_megalodon.txt -d sample_methyl.db -p 0.85

This command creates an SQLite methylation database named ``sample_methyl.db`` from Megalodon per-read methylation output,  
using a probability threshold of 0.85.




db-custom
----------

The ``db-custom`` subcommand provides a flexible interface for converting arbitrary modified basecall tables into methylartist-compatible database files.  
It is particularly useful when per-read methylation data can be exported in tabular form.

This subcommand can ingest any table that includes at least the following columns:

  - Read name  
  - Chromosome  
  - Genomic position  
  - Strand  
  - Probability that the target base is modified on a scale of 0-1.  

If the modified base is **not** specified by a column (``--modbasecol``), it can be defined manually using ``--modbase``.  
The probability column is expected to represent the **raw probability** (between 0 and 1) that a base is modified (i.e. ``1 - p(canonical)``).  
Alternative scoring schemes may be used by providing ``--canprob`` and ``--mincanprob`` to indicate a column containing canonical base scores  
and a probability cutoff for classifying canonical bases.                          

**Options:**

.. option:: -m, --methdata

   Input table containing per-read methylation data.

.. option:: --header

   Specify that the input table includes a header row.

.. option:: --delimiter

   Column delimiter character. Default is whitespace (tab or space).

.. option:: --readname

   Column number containing read names.

.. option:: --chrom

   Column number for chromosome names.

.. option:: --pos

   Column number for genomic position (0-based).

.. option:: --strand

   Column number for strand information.

.. option:: --modprob

   Column number for the probability of a base being modified.

.. option:: --canprob

   Column number for the probability of a base being canonical. If not provided, p(canonical) is assumed to be 1 - p(modified).

.. option:: --modbasecol

   Column number specifying the modified base or motif name. Can be overridden with --modbase.

.. option:: --modbase

   Manually specify the modified base or motif name. Overrides --modbasecol.

.. option:: -d, --db

   Name of the output SQLite database.

.. option:: --minmodprob

   Probability threshold for calling a base as modified. Default is 0.8.

.. option:: --mincanprob

   Probability threshold for calling a base as canonical. Default is 0.8.

.. option:: -a, --append

   Append results to an existing database.

.. option:: --motifsize

   Motif size for the modified base. Default is 2 for the common CG context.
   Set to 1 for single-base modifications such as 6mA.


**Example:**

.. code-block:: console

   methylartist db-custom -m custom_methylation.tsv --header --delimiter ',' \
       --readname 1 --chrom 2 --pos 3 --strand 4 --modprob 5 \
       --modbase CpG -d custom_methyl.db

This command loads a comma-delimited table of per-read methylation calls into an SQLite database named ``custom_methyl.db``.




db-guppy
---------

The ``db-guppy`` subcommand loads Guppy basecalled FAST5 files into a ``methylartist`` SQLite database. It is used when basecalling has been performed with Guppy and the resulting FAST5 files include modified base annotations. 

**Options:**

.. option:: -s, --samplename

   Name of the sample. Used to label entries in the database.

.. option:: -f, --fast5

   Input FAST5 file containing basecalled reads with modified base information.

.. option:: -p, --procs

   Number of processes to use for multiprocessing.

.. option:: -m, --motif

   Modified base motif (e.g. ``G[A]TC`` or ``[C]G``).

.. option:: -n, --modname

   Modified base name as defined in the Guppy FAST5 modified base alphabet (e.g. ``5mC`` or ``6mA``).

.. option:: -b, --bam

   BAM file containing alignments of reads from the FAST5 input.

.. option:: -r, --ref

   Indexed reference genome FASTA file. 

.. option:: --minprob

   Probability threshold for calling a base as modified or unmodified. Default is 0.8.

.. option:: -a, --append

   Append results to an existing database.

.. option:: --include_unmatched

   Include sites where the read base does not match the reference genome base.

.. option:: --motifsize

   Motif size for the modified base.  

.. option:: --force

   Force overwrite of existing database files if present.


**Example:**

.. code-block:: console

   methylartist db-guppy -s sample1 -f reads.fast5 -b alignments.bam -r ref.fa -m [C]G -n 5mC 

This command loads modified base calls from ``reads.fast5`` into an SQLite database,  
using ``ref.fa`` as the reference genome and labelling entries with the sample name ``sample1``.



db-sub
--------

The ``db-sub`` subcommand generates a methylation database from a **C/T-converted BAM file** that includes an ``MD`` tag.  
This command is useful when methylation status has been encoded as cytosine/thymine substitutions in the alignment file.  
It parses the BAM file, extracts modified base information using the ``MD`` tag, and stores it in an SQLite database compatible with other methylartist tools.

**Options:**

.. option:: -b, --bam

   Input BAM file containing alignments with an MD tag.

.. option:: -d, --db

   Name of the output SQLite database file. The ``.db`` extension will be added automatically if omitted.

.. option:: --append

   Append data to an existing database instead of creating a new one.



**Example:**

.. code-block:: console

   methylartist db-sub -b sample_CT.bam -d sample_CT.db

This command creates a methylation database named ``sample_CT.db`` from a C/T-converted BAM file containing an ``MD`` tag.




segmeth
-------------

``segmeth`` creates aggregate methylation or demethylation call counts over intervals. Its output is a required input for the ``segplot`` command. ``segmeth`` requires a whitespace-delimited list of segments in a BED3+2 format with chromosome, start, end, label, and strand columns. One or more BAM files may be supplied, where multiple BAM files are separated by a comma. Alternatively, the BAM files may be placed in a whitespace-delimited file with 2 fields: BAM file and .db file, which is then provided with the option ``-d``.
The bed file can be generated via ``bedtools interval``, or it can be ENCODE functional elements: promoters, enhancers, etc. It can also be used to generate genome-wide methylation statistics aggregated over a defined genomic window with

.. code-block:: console

   bedtools makewindows -g hg38.genome -w 10000 > hg38.10kbp_windows.bed

where hg38.genome is a space-delimited file with chromosome lengths that can be obtained from a fasta index file with

.. code-block:: console

   cat /path/to/.fasta.fai | awk '$2>10e6 {print $1"\t"$2}' > hg38.genome 

The output is a table showing the number of methylated, ambiguous, and unmethylated CpGs per region per sample, as well as the fraction of methylation (methylated_sites/methylated_sites + unmethylated_sites). You will find an example output table `here <https://github.com/Chisomgold/methylartist/tree/main?tab=readme-ov-file#segmeth>`__.

**Options:**

.. option:: -d, --data

   text file with .bam filename and corresponding methylation database per line (whitespace-delimited)

.. option:: -b, --bams

   one or more .bams (or bedMethyl input) with Mm and Ml tags for modification calls 

.. option:: -i, --intervals

   a .bed file

.. option:: -p, --procs

   processors for multiprocessing. Recommend for speed, if available.

.. option:: -q, --min_mapq

   minimum mapping quality (mapq) in bam file. The default is 10

.. option:: -o, --outfile

   output file name (default is generated from input)

.. option:: -m, --mods

   modifications of interest (m for 5mC); comma-delimited for >1 (default to all available modifications)

.. option:: --bedmethyl

   input is (tabix-indexed) bedMethyl

.. option:: --meth_thresh

   modified base threshold (default=0.8)

.. option:: --can_thresh

   canonical base threshold (default=0.8)

.. option:: --ctbam

   specify which .bam(s) are C/T substitution data (can be comma-delimited)

.. option:: --ref

   indexed reference genome .fa (build .fai index with samtools faidx) (required for mod bams)

.. option:: --motif

   expected modification motif (e.g. CG for 5mCpG required for mod bams)

.. option:: --max_read_density

   filter reads with call density greater >= value, can be helpful in footprinting assays (default=None)

.. option:: --excl_ambig

   do not consider reads that align entirely within segment

.. option:: --spanning_only

   only consider reads that span segment

.. option:: --primary_only

   ignore non-primary alignments

.. option:: --lowmeth_thresh

   threshold for low-methylated read count column (default = 0.05)

.. option:: --highmeth_thresh

   threshold for high-methylated read count column (default = 0.95)

.. option:: --dmr_minreads

   minimum reads per group for DMR prediction (default=8)

.. option:: --dmr_minratio

   minimum reads ratio for DMR prediction (default=0.3)

.. option:: --dmr_maxoverlap

   max overlap between distributions (default = 0.0)

.. option:: --dmr_mindiff

   minimum difference between means (default = 0.4)

.. option:: --dmr_minmotifs

   minimum motif count for DMR prediction (default = 20)

.. option:: --phased

   .bam file should be phased (e.g. with whatshap) currently only considers two phases (diploid)


**Example usage:**

.. code-block:: console

   methylartist segmeth -b first.bam,../second.bam -i hg38.10kbp_windows.bed --ref /path/to/.fasta --motif CG -p 32


segplot
------------
Use ``segplot`` to generate strip (the default), violin, or ridge plots from the segmeth output table. The category option allows for naming the elements present in the bed file used for ``segmeth``; these elements can be promoters, enhancers, histone markers, etc. You will find example plots `here <https://github.com/Chisomgold/methylartist/tree/main?tab=readme-ov-file#segplot>`_.


**Options:**

.. option:: -s, --segmeth

   output from ``methylartist segmeth``
.. option:: -m, --samples

   samples, comma-delimited
.. option:: -d, --mods

   modifications, comma-delimited
.. option:: -c, --categories

   categories, comma-delimited, need to match seg_name column from input
.. option:: -v, --violin

   generate a violin plot
.. option:: -g, --ridge

   generate a ridge plot
.. option:: -n, --mincalls

   minimum number of calls to include site (methylated + unmethylated) (default=10)
.. option:: -r, --minreads

   minimum reads in interval (default = 1)
.. option:: -q, --min_mapq

   minimum mapping quality (mapq), default = 10
.. option:: -a, --group_by_annotation

   group plots by annotation rather than by sample
.. option:: -o, --outfile

   output file name 
.. option:: --metadata

   sample metadata (tab-delimited with header, sample name as first column)
.. option:: --usemeta

   metadata to append to annotation (comma-delimited)
.. option:: --width

   figure width (default = automatic)
.. option:: --height

   figure height (default = 6)
.. option:: --pointsize

   point size for scatterplot (default = 1)
.. option:: --min,(default = -0.15)
.. option:: --max, (default = 1.15)
.. option:: --ylabel

   set label for y-axis (default: pct methylation)
.. option:: --tiltlabel

   set a title for the image
.. option:: --vertlabel
.. option:: --palette

   colour palette (default = "tab10"), see https://seaborn.pydata.org/tutorial/color_palettes.html for other options, e.g., magma
.. option:: --ridge_alpha

   transparency for ridge plot fills (default = 1.0)
.. option:: --ridge_spacing

   ridge plot spacing (generally negative, default = -0.25)
.. option:: --ridge_smoothing
   
   smoothing parameter for ridge plot, bigger is smoother (default=0.5)
.. option:: --svg

   outputs an editable vector file for Illustrator or Inkscape; useful for customised labelling for publication. Heavy files may take a while to load in Illustrator.


**Example usage:**

.. code-block:: console

   methylartist segplot -s segmeth.file.tsv -c prom,enhp -v --palette magma

This generates a violin plot showing methylation levels of promoters and enhancers in all samples present in the ``segmeth`` file.



locus
--------
Takes at least one BAM file or data file with a user-provided interval in the form of a genomic address, as chromosome:start-end. One or more regions can be highlighted using the highlight option. To generate a fully annotated image, apply gtf, fasta, and --labelgenes. Other key functionalities include variant-specific plotting and phase-based plotting (requires haplotagged BAM files, e.g., produced with whatshap).

The annotated output image and its sections from top to bottom are described `here <https://github.com/Chisomgold/methylartist/tree/main?tab=readme-ov-file#locus>`__.

**Options:**

.. option:: -d, --data

   text file with .bam filename and corresponding methylation database per line (whitespace-delimited)
.. option:: -b, --bams

   one or more .bams with MM and ML tags for modification calls
.. option:: -i, --interval

   the genomic interval to be displayed with format chrom:start-end, e.g., “chr7:1072064-1101459”. Can also be specified using flanks around the highlight region; e.g., +-10 kb means start = highlight_start - 10000 and end = highlight_end + 10000 (see -l, --highlight).
.. option:: -l, --highlight

   region(s) of interest to highlight. Format is start-end; chromosome value (if provided) is ignored and assumed to match the interval. Comma-delimited for multiple highlights.
.. option:: -g, --gtf

   genes or intervals to display in GTF format
.. option:: -C, --csi

   CSI index for the GTF (optional)
.. option:: -r, --ref

   reference genome FASTA (.fa)
.. option:: --genes

   plot and label only named genes of interest (comma-delimited)
.. option:: --labelgenes

   plot and label all genes present in the interval (requires GTF)


**Options for modifications**

.. option:: -m, --mods

   modifications of interest (e.g., m for 5mC), comma-delimited for multiple (defaults to all available modifications)
.. option:: -n, --motif

   expected modification motif (e.g., CG for 5mCpG; required for mod BAMs)
.. option:: --meth_thresh

   modified base threshold (default = 0.8)
.. option:: --can_thresh

   canonical base threshold (default = 0.8)
.. option:: --mincalls

   drop modspace positions if call count (meth + unmeth) < value (default = 0)
.. option:: --max_read_density

   filter reads with call density >= value (useful in footprinting assays; default = None)
.. option:: --excl_ambig
.. option:: --primary_only

   ignore non-primary alignments
.. option:: --unambig_highlight


**Options for coverage & masking**

.. option:: -c, --plot_coverage

   plot coverage from BAM(s) (comma-delimited list allowed)
.. option:: --logcover

   apply log2(count + 1) to coverage data (used with --plot_coverage)
.. option:: --coverprocs

   number of processes for coverage function (default = 1)
.. option:: --maskcutoff

   read-count masking cutoff (default = 1)
.. option:: --maxmaskedfrac

   skip smoothed plot if masked fraction exceeds this value (default = 1.0)
.. option:: --nomask

   skip drawing segment masks
.. option:: --slidingwindowsize

   size of initial sliding window for coverage check (default = 2)


**Options for variant, phasing & haplotypes**

.. option:: --variants

   variants to highlight (bgzipped/tabix-indexed VCF)
.. option:: --splitvar

   split reads on variant with given ID (uses ID field from VCF)
.. option:: --variantpalette

   colour palette for variant ticks (default = Set1)
.. option:: --variantsize

   size of variant tick marks (default = 6)
.. option:: --phased

   split samples into phases
.. option:: --phasediff

   add absolute difference between phases as an output track
.. option:: --ignore_ps

   ignore PS tags when plotting phased data (HP only)
.. option:: --color_by_hp

   colour samples by HP tag (requires --phased)
.. option:: --color_by_phase

   colour by phase (requires --phased)
.. option:: --include_unphased

   include an “unphased” category if plotting phased data
.. option:: --phase_labels

   custom phase labels. Format: HP:Label, comma-delimited (e.g., 1:Father,2:Mother)
.. option:: --ctbam

   specify which BAM(s) are C/T substitution data (comma-delimited)

 
**Options for smoothing**


.. option:: -s, --smoothwindowsize

   smoothing window size (default = auto)
.. option:: -t, --slidingwindowstep

   step size for initial sliding window (default = 1)
.. option:: --smoothfunc

   smoothing function: flat, hanning, hamming, bartlett, or blackman (default = hanning)
.. option:: --smoothalpha

   alpha for smoothed plot (default = 1.0)
.. option:: --smoothlinewidth

   line width for smoothed plot (default = 4.0)
.. option:: --shuffle

   shuffle sample order for smoothed plot (reduces order bias)
.. option:: --smoothed_csv

   output smoothed values to CSV (filename or "auto")


**Additional annotations**

.. option:: --bed

   BED file with additional annotations (BED3+3: chrom, start, end, label, strand, color)
.. option:: --hidebedlabel

   hide labels from BED track
.. option:: --highlight_bed

   BED3+1 with optional colour, used for multiple highlight regions


**Image formatting options:**

.. option:: --skip_align_plot

   omit read alignments and base modification glyphs
.. option:: --skip_raw_plot

   omit raw signal
.. option:: --samplebox

   draw sample box with labels next to alignments
.. option:: --allreads

   show all alignments (secondary/supplementary hidden by default)
.. option:: --readmask

   mask reads from interval(s); format start-end or chrom:start-end (chrom ignored); comma-delimited allowed
.. option:: --readmarker

   marker for (un)methylated glyphs (matplotlib markers; default = o)
.. option:: --markeralpha

   alpha for methylation markers (default = 1.0)
.. option:: --readmarkersize

   methylation marker size (default = 2.0)
.. option:: --readlinewidth

   width of alignment lines (default = 1.0)
.. option:: --readlinealpha

   transparency of alignment lines (default = 0.4)
.. option:: --readopenmarkeredgecolor

   edge colour for open (unmethylated) markers (default = sample colour)
.. option:: --exonheight

   exon height (default = 0.8)
.. option:: --show_transcripts

   plot all transcripts using transcript_id/transcript_name

.. option:: -o, --outfile

   output file name (default is autogenerated)
.. option:: --ymin

   y-axis minimum for smoothed plot (default = -0.05)
.. option:: --ymax

   y-axis maximum for smoothed plot (default = 1.05)
.. option:: --cover_ymin

   y-axis minimum for coverage plot (default = 0)
.. option:: --nticks

   number of x-axis ticks (default = 10)
.. option:: --statname

   label for raw statistics plot
.. option:: --width

   image width in inches (default = 16)
.. option:: --height

   image height in inches (default = 8)
.. option:: --notext

   remove all text from the figure
.. option:: --svg

   create vector image as output

.. option:: --colormap

   map annotations to colours; can be a file or "auto"
.. option:: --samplepalette

   palette for samples (default = "tab10")
.. option:: --coverpalette

   palette for coverage plot (default = "mako")
.. option:: --highlightpalette

   palette for highlights (default = "Blues")
.. option:: --genepalette

   palette for genes (default = "viridis")
.. option:: --highlight_alpha

   alpha value for highlight shading (default = 0.25)
.. option:: --motifsize

   motif glyph size (default = 2; set to 1 for 6mA)


**Example usage:**

.. code-block:: console

   methylartist locus -b bam1.bam,bam2.bam -i chr7:1072339-1107779 -l 1100000-1101234 -g /path/to/ref/gtf/Homo_sapiens.GRCh38.97.chr.sorted.gtf.gz --motif CG --ref /path/to/ref/hg38/Homo_sapiens_assembly38.fasta --labelgenes --width 18

This generates an annotated methylation plot of two genome samples covering a region on chromosome 7 with one highlight. The output image has a width of 18.



region
-------
Similar to ``locus`` but intended for larger regions (e.g., whole chromosomes).  
Intervals such as ``chr1:1-chr-length`` can be plotted. Chromosome lengths can be obtained from the ``.fai`` index created by ``samtools faidx``.  
Raw alignment plots are automatically omitted for intervals exceeding 5 Mb to avoid excessive rendering time; however, users can force their display if required.

**Options:**

.. option:: -i, --interval
.. option:: -d, --data

   text file with .bam filename and corresponding methylation database per line (whitespace-delimited)
.. option:: -b, --bams

   one or more BAMs with MM and ML tags for modification calls
.. option:: -l, --highlight

   highlight region(s), format start-end (chrom ignored if provided), comma-delimited for multiple regions
.. option:: -n, --motif

   normalise window sizes to motif occurrence
.. option:: -r, --ref

   reference genome FASTA (required if using motif normalisation)
.. option:: -m, --mods

   modifications to consider (comma-delimited; default = all)
.. option:: --meth_thresh

   modified base threshold (default = 0.8)
.. option:: --can_thresh

   canonical base threshold (default = 0.8)
.. option:: --ctbam

   specify which BAMs are C/T-substitution data (comma-delimited)
.. option:: --bedmethyl

   input is tabix-indexed bedMethyl
.. option:: --max_read_density

   filter reads with call density ≥ value (default = None)

.. option:: -c, --plot_coverage

   plot coverage from BAMs (comma-delimited allowed)
.. option:: -w, --windows

   number of windows to use (default = auto)
.. option:: --maxuncovered

   maximum percentage of uncovered windows tolerated (default = 50)
.. option:: --min_window_calls

   minimum reads per window (default = 1)
.. option:: -q, --min_mapq

   minimum mapping quality (default = 10)
.. option:: -p, --procs

   multiprocessing
.. option:: --logcover

   apply log2(count+1) to coverage values

.. option:: --readmask

   mask reads in interval(s); format start-end or chrom:start-end, comma-delimited
.. option:: --allreads

   show all alignments (including secondary/supplementary)
.. option:: --primary_only

   ignore non-primary alignments


**Smoothing options**

.. option:: -s, --smoothwindowsize
.. option:: --smoothfunc

   one of: flat, hanning, hamming, bartlett, blackman (default = hanning)
.. option:: --smoothalpha

   transparency for smoothed plot (default = 1.0)
.. option:: --smoothlinewidth

   smoothed plot line width (default = 4.0)
.. option:: --shuffle

   shuffle sample order in smoothed plot
.. option:: --segment_csv

   export smoothed window/segment values to CSV


**Annotation options**

.. option:: -g, --gtf

   GTF file for gene annotation
.. option:: -C, --csi

   CSI index for GTF (optional)
.. option:: --genes

   plot only named genes (comma-delimited)
.. option:: --labelgenes

   label gene names
.. option:: --show_transcripts

   plot all transcripts (uses transcript_id / transcript_name)
.. option:: --exonheight

   exon height (default = 0.8)
.. option:: --gene_track_height

   maximum number of gene track layers
.. option:: --bed

   BED file for extra annotations (BED3+3)
.. option:: --hidebedlabel

   hide labels from BED track
.. option:: --highlight_bed

   BED3+1 with optional colour for highlight regions


**DMR-related analysis**

.. option:: --dmr_minreads

   minimum reads per group (default = 8)
.. option:: --dmr_minratio

   minimum read ratio (default = 0.3)
.. option:: --dmr_maxoverlap

   max overlap between distributions (default = 0.0)
.. option:: --dmr_mindiff

   minimum difference between mean values (default = 0.4)
.. option:: --dmr_minmotifs

   minimum motif count (default = 20)
.. option:: --write_dmrs

   output differentially methylated windows

**Phasing options**

.. option:: --phased

   currently assumes diploid (two phases)
.. option:: --color_by_hp

   colour samples by HP value (requires --phased)

**Image formatting options:**

.. option:: --skip_align_plot

   blank alignment panel (useful for long regions)
.. option:: --force_align_plot

   force alignment plots even for >5 Mbp regions
.. option:: --modspace

   spacing between links in top panel (default = auto)
.. option:: --motifsize

   motif glyph size (default = 2; set to 1 for 6mA)
.. option:: --panelratios

   four (or five with coverage) comma-separated integers defining panel height ratios (default: 1,5,3,3)

.. option:: --highlight_alpha

   highlight transparency (default = 0.25)
.. option:: --highlight_centerline

   draw highlight as a center line with specified width
.. option:: --samplepalette

   palette for samples (default = "tab10")
.. option:: --coverpalette

   palette for coverage plot (default = "mako")
.. option:: --highlightpalette

   palette for highlights (default = "Blues")
.. option:: --genepalette

   palette for gene annotation (default = "viridis")
.. option:: --colormap

   annotation-to-colour mapping file or "auto"
.. option:: --ymin
.. option:: --ymax
.. option:: --cover_ymin

   y-axis minimum for coverage (default = 0)
.. option:: --nticks

   tick count (default = 10)
.. option:: --scale_fullwidth

   scale plot width relative to value (e.g., chromosome length)
.. option:: --width

   image width in inches (default = 16)
.. option:: --height

   image height in inches (default = 8)
.. option:: --svg
.. option:: -o, --outfile

   output filename (default = autogenerated)


**Example usage:**

.. code-block:: console

   methylartist region -i chr22:1-50818468 -b first.bam,second.bam -p 32 -g /path/to/.gtf.gz -n CG -r /path/to/.fasta

This generates a methylation plot across the entire chromosome 22 for two samples. Alignment plots are disabled automatically because the region exceeds 5 Mb.



composite
--------------
Generates composite methylation plots, where multiple per-element profiles are aligned to and plotted against a reference sequence.  
An example input file and output image are available `here <https://github.com/Chisomgold/methylartist/tree/main?tab=readme-ov-file#composite>`__.

**Options:**

.. option:: -d, --data

   text file with .bam filename and corresponding methylation database per line (whitespace-delimited)
.. option:: -b, --bams

   one or more BAMs with MM and ML tags for modification calls
.. option:: -s, --segdata

   BED3+1 file with chrom, start, end, strand
.. option:: -r, --ref

   reference genome FASTA
.. option:: -t, --teref

   TE reference FASTA
.. option:: -p, --procs

   number of processes to use

.. option:: --mod

   modification to plot (default inferred from sample name)
.. option:: --motif

   motif to highlight (default = CG)
.. option:: --meth_thresh

   modified base threshold (default = 0.8)
.. option:: --can_thresh

   canonical base threshold (default = 0.8)
.. option:: -q, --min_mapq

   minimum mapping quality (default = 10)
.. option:: -l, --lenfrac

   minimum fraction of TE length that must align (default = 0.95)
.. option:: --mincalls

   minimum call count to include an element (default = 100)
.. option:: --minelts

   minimum number of elements required in output (default = 1)
.. option:: --maxelts

   maximum number of output elements; if exceeded, random sampling is used (default = 200)
.. option:: --max_read_density

   filter reads with call density ≥ value (default = None)
.. option:: --excl_ambig
.. option:: --primary_only

   ignore non-primary alignments
.. option:: --phased

   assumes two phases (diploid)
.. option:: --color_by_phase

   colour samples by HP value (requires --phased)

.. option:: --slidingwindowsize

   size of sliding window for methylation fraction (default = 10)
.. option:: --slidingwindowstep

   step size for sliding window (default = 1)
.. option:: --smoothwindowsize

   smoothing window size (default = 8)
.. option:: --smoothfunc

   smoothing function: flat, hanning, hamming, bartlett, blackman (default = hanning)

.. option:: --start

   start plotting at this base (default = None)
.. option:: --end

   end plotting at this base (default = None)
.. option:: --blocks

   blocks to highlight (text file: start, end, name, hex colour)
.. option:: --meanplot_ylabel

   set y-axis label for mean plot
.. option:: --meanplot_cutoff

   override site coverage cutoff for mean plot (automatic value shown in output)

.. option:: -c, --palette

   palette for samples (default = "tab10")
.. option:: -a, --alpha

   transparency (default = 0.3)
.. option:: -w, --linewidth

   line width (default = 1)
.. option:: --ymin
.. option:: --ymax
.. option:: -o, --outfile

   output filename (default autogenerated)
.. option:: --output_table

   output per-site data as a TSV table
.. option:: --svg

   produce an SVG (vector) figure


wgmeth
--------

Takes a BAM file, a FASTA index file, the FASTA file, the motif (default CG), modification type, number of processors, and optional DSS output mode. Produces genome-wide methylation statistics for all CpGs (or motif sites) covered by at least one call.

By default, output is in bedMethyl format (0-based). Alternatively, output can be a 1-based table compatible with DSS for DMR calling on nanopore data.

**Options:**

.. option:: -b, --bam

   BAM used for methylation calling
.. option:: -d, --methdb

   methylation database
.. option:: -r, --ref

   reference genome FASTA (required for mod BAMs; build index with *samtools faidx*)
.. option:: -f, --fai

   FASTA index (.fai). Default = <ref>.fai. Required for .db files.
.. option:: --motif

   motif to query (e.g., CG for 5mCpG; required for mod BAMs)
.. option:: -m, --mod

   modification code to output (names vary; see output for hints)
.. option:: -p, --procs

   number of processes
.. option:: -s, --binsize

   bin size for parallelisation (default = 1,000,000)
.. option:: -c, --chrom

   limit analysis to one chromosome
.. option:: -l, --minlen

   minimum chromosome length (default = 0)
.. option:: -q, --min_mapq

   minimum mapping quality (default = 10)
.. option:: --meth_thresh

   modified base threshold (default = 0.8)
.. option:: --can_thresh

   canonical base threshold (default = 0.8)
.. option:: --max_read_density

   filter reads with call density ≥ value (default = None)
.. option:: --ctbam

   specify which BAMs are C/T substitution data (comma-delimited)
.. option:: --dss

   output in DSS format (1-based) instead of bedMethyl (0-based)
.. option:: --phased

   split output into phases (currently 1 and 2)
.. option:: --primary_only

   ignore non-primary alignments
.. option:: -o, --outfile

   output filename (default generated from input)


**Example usage:**

.. code-block:: console

   methylartist wgmeth -b file.bam -f /path/to/fasta.fai -r /path/to/fasta \
       --motif CG --dss --mod m -p 32

This generates a table with the following information for each CG position in the query genome sample: chromosome number, genomic coordinate, total number of reads, and number of reads showing methylation.



adjustcutoffs
-------------------
This is used to update a database file (in-place) with user-defined cut-offs for methylated (1) state and unmethylated (0) state. 

**Options:**

.. option:: -d, --db

   methylartist database
.. option:: --mod

   modification to plot (will list for user if incorrect)
.. option:: -m, --methylated

   mark as methylated above cutoff value
.. option:: -u, --unmethylated

   mark as unmethylated below cutoff value


scoredist
---------
Outputs the distribution of modified base statistics using current cutoffs for one or more methylartist databases.

**Options:**

.. option:: -d, --db

   methylartist database(s), comma-delimited
.. option:: -b, --bam

   one or more BAM files with MM/ML tags for modification calls
.. option:: -m, --mod

   modification to plot (will list available options if incorrect)
.. option:: --motif

   modified motif to highlight (e.g., CG)
.. option:: -r, --ref

   reference genome FASTA (samtools faidx indexed)
.. option:: -n, --n

   sample size (default = 1,000,000)
.. option:: --xmin

   minimum x-axis value
.. option:: --xmax

   maximum x-axis value
.. option:: --lw

   line width (default = 2)
.. option:: --palette

   colour palette for phases (default = "tab10")
.. option:: -o, --outfile

   output file name (default generated from input)
.. option:: --svg

   output SVG format


**Example:**

.. code-block:: console

   methylartist scoredist -d MCF7_ATCC.megalodon.db,MCF7_ECACC.megalodon.db -m m

