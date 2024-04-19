##' polyG
##' 
##' Main script to process polyG data, from raw data input to angular distance data.
##' 
##' @export 
polyG <- function(input_dir, results_dir, max_replicate_th=0.11, sample_name_change=T, bootstrap=T, bscut=50, bsreps=1000, seed=NA) {  

    # load the data_info file for this patient
    data_info <- read.table(file.path(input_dir,"data_info.txt"),sep="\t",header=TRUE,stringsAsFactors=FALSE)

    # process to run for each patient
    for (prow in 1:nrow(data_info)) {
        if(!is.na(seed)) set.seed(seed) ## reset seed for each patient to enable rerunning for single patients at a time

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
        new_sample_names <- read.table(file.path(input_dir,"annotation_files",paste0(subject_name,"_annotationFile.txt")),sep="\t",header=TRUE,stringsAsFactors=FALSE)

        # create output directory
        allout_dir <- file.path(results_dir,"results",paste("sample_exclusion",sample_exclusion_th,"rep_cut",max_replicate_th,sep="_"))
        if (!dir.exists(allout_dir)) dir.create(allout_dir,recursive=TRUE,showWarnings=TRUE )

        # select representative replicates and create distance matrices
        subject(input_dir, allout_dir, subject_name, data_path, max_replicate_th, sample_exclusion_th, sel_normal_sample, all_normal_samples, new_sample_names) 

        ## generate angular distance matrices, heatmaps, trees, and optionally bootstrapped trees
        angular_distance(input_dir, allout_dir, subject_name, sel_normal_sample, all_normal_samples, new_sample_names, sample_name_change, bootstrap, bscut, bsreps)
    }
}


