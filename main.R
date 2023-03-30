### This script is used to analyze poly-G data for all patients.
### Requires: 
### 2. for each patient, all marker files in a directory "input_dir/data"
### 3. for each patient, a column of all sample names in a file named "patient ID" + "_SampleNames.txt" in directory "input_dir/sample_names"
### 4. A file in directory "input_dir/colorfiles" for tree tip colors with filename of the form of "patient" + "_colorfile.txt" (e.g., E18_colorfile.txt)
### 5. A file in directory "input_dir/AnnotationFiles" with filename of the form of "patient" + "_annotationFile.txt") (e.g., E1_annotationFile.txt) 
###    The annotation files should contain two columns with column names. The first column should be the original sample names 
###      and the second column the final sample names (which can be the same as the original sample names.)
### 6. data_info.txt which contains:
###    (a) one patient per row
###    (b) first column: data_path for a patient, e.g., "E1-Data"
###    (c) second column: sel_normal_sample (the normal samples used as reference, e.g., "E1N2")
###    (d) 3rd to (nnormals+1)th columns: list other normal samples used in create_marker_files step, one per column; 
###         nnormals = maximum number of normal samples among all the patients
###    (e) additonal columns: list all other samples to exclude for KLM calculation, 
###            one sample per column (e.g.,"E1H1b" - keep only one sample per lesion)
### 7. adjust all relevant parameters (including a "nnormals") above the line with ########


## install any missing R CRAN packages:
required.R.packages <- c('here','gplots','RColorBrewer','colorspace')
new.R.packages <- required.R.packages[!(required.R.packages %in% installed.packages()[,"Package"])]
if(length(new.R.packages)>0) install.packages(new.R.packages)

## install any missing R Bioconductor packages:
required.BC.packages <- c('phyloseq','dendextend','ape','phangorn','adephylo')
new.BC.packages <- required.BC.packages[!(required.BC.packages %in% installed.packages()[,"Package"])]
if(length(new.BC.packages)>0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(new.BC.packages)
}

## load required packages
required.packages <- c(required.R.packages, required.BC.packages)
for(p in required.packages) library(package=p, character.only=T)

## load functions
source(here('polyG_pipeline/func.R'))
source(here('polyG_pipeline/subject.R'))
source(here('polyG_pipeline/angular_distance.R'))


# the folder with the input data for this patient cohort
input_dir <- here("original_data/polyG_lung")

# where the output "results" will go
results_dir <- here('processed_data/polyG_lung')

# the maximum JS distance allowed between replicates
max_replicate_th <- 0.11

# plotting heatmaps and trees with different sample names from those in the marker files? 
sample_name_change <- TRUE

# bootstrap options
bootstrap <- F
bscut <- 50 ## not currrently implemented

# load the data_info file for this patient
data_info <- read.table(file.path(input_dir,"data_info.txt"),sep="\t",header=TRUE,stringsAsFactors=FALSE)


#####################
# for each patient
for (prow in 1:nrow(data_info)) {
    set.seed(42) ## reset seed for each patient to enable rerunning for single patients at a time
 
    subject_name <- gsub('-Data','',data_info$data_path[prow])
    message('Running for ',subject_name)

    ## get path, exclusion threshold and selected normal for this patient
    data_path <- data_info$data_path[prow]
    sample_exclusion_th <- data_info$sample_exclusion_th[prow]
    sel_normal_sample <- data_info$sel_normal_sample[prow]

    ## get all the normal samples for this patient
    all_normal_samples <- as.character(data_info[prow,grep('normal_sample',names(data_info))])
    all_normal_samples <- all_normal_samples[!is.na(all_normal_samples)]

    ## get all excluded samples for this patient
    all_excluded_samples <- as.character(data_info[prow,grep('samples_to_exclude',names(data_info))])
    all_excluded_samples <- all_excluded_samples[!is.na(all_excluded_samples)]
  
    ## map old/new names 
    new_sample_names <- read.table(file.path(input_dir,"AnnotationFiles",paste0(subject_name,"_annotationFile.txt")),sep="\t",header=TRUE,stringsAsFactors=FALSE)

    # create output directory
    allout_dir <- file.path(results_dir,"results",paste("sample_exclusion",sample_exclusion_th,"rep_cut",max_replicate_th,sep="_"))
    if (!dir.exists(allout_dir)) dir.create(allout_dir,recursive=TRUE,showWarnings=TRUE )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # select representative replicates and create distance matrices
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subject(input_dir, allout_dir, subject_name, data_path, sample_exclusion_th, sel_normal_sample, all_normal_samples, new_sample_names) 
    
    ## get angular distance trees
    angular_distance(input_dir, allout_dir, subject_name, sel_normal_sample, all_normal_samples, new_sample_names, bootstrap, bscut)
}




