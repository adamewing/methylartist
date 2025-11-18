.. methylartist documentation master file, created by
   sphinx-quickstart on Mon Nov 17 13:35:36 2025.

Methylartist:  A suite of tools for visualising nanopore-derived modified bases data
====================================================================================

Methylartist is a Python-based toolkit for visualising methylation data.
It provides utilities to create plots, including violin and ridge plots, from nanopore-derived sequence alignment files (bam) for either the whole genome or regions of interest. In addition, it supports the generation of DSS-compatible file types, modification database files (.db), and others.

Key Features
------------

- Loci-specific methylation plotting
- Simple command line interface
- Whole chromosome methylation visualisation

Citation
--------

If you use methylartist in your research, please cite:

.. code-block:: console

   @article{10.1093/bioinformatics/btac292,
    author = {Cheetham, Seth W and Kindlova, Michaela and Ewing, Adam D},
    title = {Methylartist: tools for visualizing modified bases from nanopore sequence data},
    journal = {Bioinformatics},
    volume = {38},
    number = {11},
    pages = {3109-3112},
    year = {2022},
    month = {04},
    abstract = {Methylartist is a consolidated suite of tools for processing, visualizing and analysing nanopore-derived modified base calls. All detectable methylation types (e.g. 5mCpG, 5hmC, 6mA) are supported, enabling integrated study of base pairs when modified naturally or as part of an experimental protocol.Methylartist is implemented in Python and is installable via PyPI and bioconda. Source code and test data are available at https://github.com/adamewing/methylartist.Supplementary data are available at Bioinformatics online.},
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btac292},
    url = {https://doi.org/10.1093/bioinformatics/btac292},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/38/11/3109/49878586/btac292.pdf},
   }


License
--------

This project is licensed under the MIT License.

Documentation by Halimat Chisom Atanda

Date: Oct 13, 2025

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   commands
