##' sel_rep_replicates_create_dist_matrix
##' @export
sel_rep_replicates_create_dist_matrix <- function(subject_name,sel_normal_sample,all_normal_samples,pdir='.',overwrite=F,sample_exclusion_th=0.15,max_replicate_th=0.11,intfraction_cutoff=0.1,abovebelow_cutoff=0.15,distcut=0.06,file_filter=0.1,rob_th=0.05,sample_name_change=T) {
    ## select representative replicates and create distance matrices for each sample from the given subject

    ## define the path to the various poly-G loci data files for this given subject
    data_path <- file.path(pdir,'data',paste0(subject_name,'-Data'))
    data_files <- list.files(data_path)

    ## define the expected output directories
    date <- format(Sys.Date(),format="%Y%m%d")
    allout_dir <- file.path(pdir,"results",paste("sample exclusion",sample_exclusion_th,"infraction",intfraction_cutoff,"abovebelow",abovebelow_cutoff,"distcut",distcut,"rep_cut",max_replicate_th,date,sep="_"))
    output_dir <- file.path(allout_dir,paste0(subject_name,"_",file_filter,"_",date,"_R"))

    ## check if expected output already exists before running. If so, and overwrite==F, throw error.
    create_dir_check_if_already_exists(allout_dir,overwrite)
    create_dir_check_if_already_exists(output_dir,overwrite)

    ## run for subject, collect output from this process to then generate the distance_distribution_heatmap
    subject_output_list <- subject(subject_name,pdir,output_dir,sample_exclusion_th,max_replicate_th,sel_normal_sample,all_normal_samples)
    marker_list <- subject_output_list$marker_list
    sample_names <- subject_output_list$sample_names   

    if(!sample_name_change) { 
        ## run with only the single set of sample names
        for (markers in marker_list) {
            ddh_fp <- get_output_filepath(subject_name,"heatmap",paste0("_",markers,".pdf"),output_dir,max_replicate_th)
            distance_distribution_heatmap(ddh_fp,subject_name,markers,max_replicate_th,sel_normal_sample,all_normal_samples,rob_th,subject_output_list)
        }

    } else {
        ## user specified to use new sample names, load the _annotationFile.txt to extract these.
        new_rownames <- sample_names ## start with default sample names
        if (sample_name_change) {
            anno_file <- file.path(pdir,"AnnotationFiles",paste0(subject_name,"_annotationFile.txt"))
            error_if_file_not_exist(anno_file)
            anno <- read.table(anno_file,sep="\t",header=T,stringsAsFactors = F)
            new_rownames <- anno[pmatch(sample_names,anno[,1]),2]
        }

        for (markers in marker_list) {
            ddh_fp_old <- get_output_filepath(subject_name,"heatmap",paste0("_",markers,"_oldnames.pdf"),output_dir,max_replicate_th)
            distance_distribution_heatmap(ddh_fp_old,subject_name,markers,max_replicate_th,sel_normal_sample,all_normal_samples,rob_th,subject_output_list)
            ddh_fp_new <- get_output_filepath(subject_name,"heatmap",paste0("_",markers,"_newnames.pdf"),output_dir,max_replicate_th)
            distance_distribution_heatmap(ddh_fp_new,subject_name,markers,max_replicate_th,sel_normal_sample,all_normal_samples,rob_th,subject_output_list,new_rownames=new_rownames)
        }
    }

    #if (length(final_marker_names[["used_markers"]]) > 0) {
    #    source("polygR/create_distance_matrix.R")
    #    source("polygR/create_clustermap.R")
    #} else {
    #    print("No usable markers available. Clustermap cannot be generated.")
    #}

}
