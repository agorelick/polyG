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


library(polyG)

# the folder with the input data for this patient cohort
input_dir <- '~/lab_repos/crc_lung_met/original_data/polyG_lung'

# where the output "results" will go
results_dir <- '~/lab_repos/crc_lung_met/processed_data/polyG_lung'

# the maximum JS distance allowed between replicates
max_replicate_th <- 0.11

# plotting heatmaps and trees with different sample names from those in the marker files? 
sample_name_change <- TRUE

# bootstrap options
bootstrap <- F
bscut <- 50 ## not currrently implemented

polyG(input_dir, results_dir, max_replicate_th, sample_name_change, bootstrap, bscut, seed=42)



