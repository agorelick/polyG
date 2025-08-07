# polyG

## Overview

This GitHub repo contains the _polyG_ R package, which implements the pipeline used to process raw polyguanine genotype data as described in _Wassenaar, Gorelick, et al. 2024_. The instructions described below can additionally be used to process new polyguanine genotype data, provided the input files is identically prepared. For analyses and data related to  _Wassenaar, Gorelick, et al. 2024_, please see the accompanying GitHub repo: https://github.com/agorelick/peritoneal_metastasis. 

## 0. Clone this GitHub repo

On Mac/Linux, enter the following in a terminal:
```bash
cd <location to download this GitHub repo>
git clone https://github.com/agorelick/polyG
cd polyG
```

## 1a. Installation (using conda)

```bash
# create a new conda environment (here called "polyG") and install R and necessary packages into it
conda env create -n polyG -f env.yml -y

# activate it
conda activate polyG

# install the polyG R package into it
R CMD INSTALL . 
```

## 1b. Installation without conda (using your default R installation)

Open R/Rstudio, install prerequsite R packages manually
```r
# Install missing R packages from CRAN:
required.R.packages <- c('here','gplots','RColorBrewer','colorspace')
new.R.packages <- required.R.packages[!(required.R.packages %in% installed.packages()[,"Package"])]
if(length(new.R.packages)>0) install.packages(new.R.packages)

# Install missing R packages from Bioconductor:
required.BC.packages <- c('phyloseq','dendextend','ape','phangorn','adephylo')
new.BC.packages <- required.BC.packages[!(required.BC.packages %in% installed.packages()[,"Package"])]
if(length(new.BC.packages)>0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(new.BC.packages)
}
```

Next, install the polyG package:

```r
install.packages('[PATH TO CLONED REPO]/polyG/pkg', type='src', repos=NULL)
```


## 2. Run the poly-G pipeline on the provided example data:

This pipeline is used to analyze poly-G data for all patients. The **following input data are required** (for examples, see the files in the "example_input/" directory): 

1. For each patient, all marker files in a directory (e.g., "example_input/data")
2. For each patient, a file with sample names (corresponding the sample names provided in the raw marker files data) to be used for each patient; samples not listed will be excluded from analysis (e.g., "example_input/C31_SampleNames.txt")
3. For each patient, a file with tip colors for each sample according to various features, such as "tissue_type" (e.g., "example_input/color_files/C31_colorfile.txt")
4. For each patient, a file mapping "raw" sample names (used in the marker files in step 1) to "finalized" sample names. Note that if new sample-names are not required, just use the same values in both columns of these files (e.g., "example_input/annotation_files/C31_annotationFile.txt")
5. **data_info.txt** - a tab-delimited table containing one patient per row, with the following columns:
    - Column 1: "data_path" - the filename for each patient's raw marker file (e.g., "C31-Data")
    - Column 2: "sel_normal_sample" - a single normal sample (using the "raw" sample names) to use as a germline reference for each patient (e.g., "C31-Y")
    - Columns 3-5: list any other normal samples for patient, used in the create_marker_files step. If any sample has more than 3 normal samples, just add more "other_normal_sample_" columns to the table following the naming convention. (e.g., "C36-1-C" for patient C36)
    - Columns 6-12: list any samples to exclude for each patient (using the "raw" sample names, e.g., "C31-M" for patient C31)

   
```r
# load the polyG R library
library(polyG)

# run poly-G processing pipeline for the small example dataset (patients C31, C36 from Naxerova et al, 2017 Science)
polyG(input_dir='example_input', results_dir='example_output', seed=42)
```
