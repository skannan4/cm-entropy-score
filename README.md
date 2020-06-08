# Transcriptomic entropy quantifies cardiomyocyte maturation at single cell level

### Introduction
Transcriptomic entropy, quantified from single cell RNA-sequencing (scRNA-seq) data, can be used to quantify the maturation status of cardiomyocytes (CMs). This may be particularly helpful for assessing the maturation status of pluripotent stem cell-derived CMs (PSC-CMs) and other engineered CM tissues. The rationale and details of our approach can be found in our preprint (add citation when ready). Here, we provide all of the relevant code and workspaces for replicating our figures in the preprint and for using entropy score to analyze your own datasets.

### Data
The following relevant files can be downloaded at our [Synapse](https://www.synapse.org/#!Synapse:syn21788425/files/):

- `clean_060720.RData`: An R workspace containing all of the counts tables for all datasets analysed in our manuscript, the computed QC metrics/entropy scores/other relevant metadata for all datasets, and all functions used in our manuscript. Please note that this file is very large and may require significant RAM to load and work off of. We recommend 32GB RAM for handling this file.

- `clean_nodatasets_060720.RData`: An R workspace containing all of the computed QC metrics/entropy scores as well as functions used in the manuscript, but not the counts tables. This file cannot be used to reproduce every figure in the manuscript. However, it may be useful for those with lower RAM availability, and can be used to run entropy score on user datasets.

- `supp_files.zip`: Zip file containing two folders (Fig02_data and Fig04_data) with some necessarily files for replicating a handful of the figures.

The following relevant files can be found here:


- `entropy_functions.R`: An R file containing code for all of the relevant functions for entropy score. Please note that these functions have been loaded into both of the above workspaces as well.

- `entropy_figures.R`: A R file containing code necessary to reproduce all figures in the manuscript. Please note that many of the figures will require `clean_060720.Rdata`.

- `helper_code.R`: An R file containing some miscellaneous helper codes that were beneficial throughout the project.

### Important objects in the workspaces
In addition to the counts tables in `clean_060720.Rdata`, the workspaces have several other objects of interest. Many are related specifically to the functions, and their use is detailed in `entropy_functions.R`. However, two dataframes may be of particular interest to users:

- `alldata`: This contains the metadata for every dataset that we analyzed. It was also used to generate the Supplementary Tables in the manuscript. Readers looking for details about datasets should look here.

- `combined_datasets`: This is the primary output dataframe of our workflow. It contains the tabulated quality control metrics, computed entropy scores, and all other relevant details for every dataset analyzed. For those looking solely for the downstream output of our workflow, this may be appropriate dataframe.

### Dependencies
Most of the libraries used in our codebase can be found from CRAN or Bioconductor. However, we additionally make use of the SingleCellNet package from the Cahan lab. Please see [their github](https://github.com/pcahan1/singleCellNet) for instructions on how to install SingleCellNet.

### How to replicate figures in our manuscript
To replicate the figures in our manuscript, please follow these steps:

1. Download `clean_060720.Rdata`, `entropy_figures.R`, and `supp_files.zip`.

2. Load `clean_060720.Rdata` into R or your IDE of choice.

3. Extract `supp_files.zip` to a folder of your choice.

4. Modify the first line of `entropy_figures.R`: `setwd("~/Documents/Research/Reference")` to set the working directory to the same working directory you chose in Step 3.

5. Load in all the required packages listed at the top of `entropy_figures.R`.

6. Run code for any figure of interest. The figures are broken into self-contained code-blocks; at the end of each code block, all temporary objects created for that figure are removed. Each code block should also have details on the figure itself that may be helpful.

Please feel free to email or raise an issue if any of the code doesn't work as claimed!

### How to use entropy score for your own datasets
Entropy score is relatively straightforward to use and should extend robustly to many datasets (with minor limitations, particularly at very low depth, noted in the manuscript). If you are interested in testing entropy score for your own datasets, please follow these steps:

1. Make sure you have SingleCellNet installed.

2. Download one of `clean_060720.Rdata` or `clean_nodatasets_060720.Rdata` and load into R - all of the necessary functions have been already loaded in these workspaces. The latter may be helpful for those with memory limitations or just looking to rapidly test entropy score on their dataset. If you would like more details about each function, please download `entropy_functions.R`.

3. Load and format your dataset appropriately. You will need, at minimum, the following: a counts table where the rows are gene names (in gene symbol format) and the columns are cellnames; a vector containing all of the timepoints of your cells (in the same order as your columns in the counts table, but doesn't need to be named). Other metadata about the study may also be helpful, but is unnecessary.

  * If your dataset has rownames in ENSEMBL format, you may use the `rename_genes()` function. For example, for mouse datasets, run `rename_genes(dataset, species = "mouse")`; for human datasets, run, `rename_genes(dataset, species = "human")`.

4. If you are 1000% percent convinced that your dataset has been appropriately filted for CMs and only high quality cells, you may simply use the `master_entropy()` function - e.g. `master_entropy(dataset)`. This will return a vector containing the entropy score for each cell in your dataset. However, we strongly recommend you use our QC approach, outlined in Step 5.

5. We provide the function `data_qc()`, which serves as a one-shot function for getting all of the QC information necessary for appropriately using entropy score. `data_qc()` has the following parameters (please skip to step six if you just want sample code on how to run):
  * `dataset`: REQUIRED. Your dataset object name as a character. Please note that, unlike the other functions, your dataset name MUST BE IN QUOTES.
  * `study`: A character giving the name of the study. (e.g. "Kannan et al.")
  * `timepoint_list`: REQUIRED. Vector listing the timepoints of all the cells in the dataset.
  * `scn_calc`: Whether to run the SingleCellNet function to determine whether your cells are CMs. Skipping will significantly speed up this step (as SingleCellNet is the most time-consuming portion of `data_qc()`); however, we strongly urge you run this step even if you are confident in the identity of your cells.
  * `species`: Can be set to "mouse" or "human".
  * `sample_type`: Optional. I use to store "in vivo", "directed differentiation", or "direct reprogramming".
  * `isolation`: Optional. I use to store the method by which sample was acquired, e.g. FACS, manual picking.
  * `sequencing`: Optional. I use to store the sequencing protocol, e.g. SCRB-seq.
  * `mapping`: Optional. I use to store the mapping method, e.g. STAR/FeatureCounts, zUMIs, kallisto.
  * `datatype`: Optional. I use to store the datatype, e.g. reads, UMIs.
  * `doi`: Optional.
  * `other_meta`: Optional. If you have another metadata field of interest (for example, atrial vs. ventricular), you can input here as a vector, much as with timepoint.

The output of this will be a dataframe that has all of the relevant information for using entropy score. In particular, the following fields will be useful:
  * `max_celltype`: We select only cells labelled as "cardiac muscle cell."
  * `good_cell`: Set to true if `top5_norm` < 1.3 and `depth_norm` > -0.5.
  * `genes`: We select cells > 1000.
  * `entropy`: This is what you are here for!

6. Here is an example snippet:
```R
 wonderful_temp = data_qc(dataset = "wonderful_data", study = "Brilliant Grad Student et al.", timepoint_list = wonderful_timepoint_list, scn_calc = TRUE, species = "mouse", sample_type = "in vivo", isolation = "LP-FACS", sequencing = "mcSCRB-seq", mapping = "zUMIs", datatype = "UMIs", doi = "doi:12345", other_meta = NA)
 ```
 
 7. Analyze your data however you choose (for example in `ggplot2`). Note - it is incredibly straightforward to analyze your data alongside all of the other datasets in our meta-analysis. Simply run `combined_datasets = rbind(combined_datasets, wonderful_temp`, where `wonderful_temp` was the output from Step 6. In fact, `combined_datasets` was actually constructed exactly this way, running `data_qc()` on every dataset we found!
 
 8. Exclaim with glee as you have now quantified your CM maturation status.
