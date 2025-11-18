.. _installation:

Installing methylartist
=======================

Methylartist can be installed through pip, conda/bioconda or the repository. 

Installing with pip or conda
-----------------------------

To install methylartist in your current conda environment, type: 

``pip install methylartist`` OR
 
``conda install bioconda::methylartist`` OR 

``conda install -c bioconda methylartist``

This ensures you automatically install all dependencies and potentially avoid OS issues.

Installation from the repository
---------------------------------

Clone the repository with 

.. code-block:: console

   git clone https://github.com/adamewing/methylartist.git
   cd methylartist
   conda create methylartist 
   conda activate methylartist
   pip install -e $PWD
   methylartist -h

This should display a help page, confirming that the package is successfully installed. In case of an update to the repo: 

.. code-block:: console

   cd path/to/methylartist
   git pull
   pip install -e $PWD

