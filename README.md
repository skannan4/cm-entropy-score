# Transcriptomic entropy quantifies cardiomyocyte maturation at single cell level

### Introduction
Transcriptomic entropy, quantified from single cell RNA-sequencing (scRNA-seq) data, can be used to quantify the maturation status of cardiomyocytes (CMs). This may be particularly helpful for assessing the maturation status of pluripotent stem cell-derived CMs (PSC-CMs) and other engineered CM tissues. The rationale and details of our approach can be found in our preprint (add citation when ready). Here, we provide all of the relevant code and workspaces for replicating our figures in the preprint and for using entropy score to analyze your own datasets.

### Data
The following relevant files can be downloaded here:

- clean_031020.RData: An R workspace containing all of the counts tables for all datasets analysed in our manuscript, the computed QC metrics/entropy scores/other relevant metadata for all datasets, and all functions used in our manuscript. Please note that this file is very large and may require significant RAM to load and work off of. We recommend 32GB RAM for handling this file.

- clean_nodatasets_031020.RData: An R workspace containing all of the computed QC metrics/entropy scores as well as functions used in the manuscript, but not the counts tables. This file cannot be used to reproduce every figure in the manuscript. However, it may be useful for those with lower RAM availability, and can be used to run entropy score on user datasets.

- entropy_functions.R: An R file containing code for all of the relevant functions for entropy score. Please note that these functions have been loaded into both of the above workspaces as well.

- entropy_figures.R: A R file containing code necessary to reproduce all figures in the manuscript. Please note that many of the figures will require clean_031020.Rdata.

- helper_code.R: An R file containing some miscellaneous helper codes that were beneficial throughout the project.
