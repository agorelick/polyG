##' create_phylogenetic_trees_with_bootstrap_values
##' @export
create_phylogenetic_trees_with_bootstrap_values <- function(subject_name,sel_normal_sample,all_normal_samples,pdir,overwrite,sample_exclusion_th,max_replicate_th,intfraction_cutoff=0.1,abovebelow_cutoff=0.15,distcut=0.06,file_filter=0.1,rob_th=0.05,sample_name_change=T,colorfile_header=F,KLM=F,bscut=50) {

    ## define the path to the various poly-G loci data files for this given subject
    data_path <- file.path(pdir,'data',paste0(subject_name,'-Data'))
    data_files <- list.files(data_path)

    ## define the expected output directories
    date <- format(Sys.Date(),format="%Y%m%d")
    allout_dir <- file.path(pdir,"results",paste("sample exclusion",sample_exclusion_th,"infraction",intfraction_cutoff,"abovebelow",abovebelow_cutoff,"distcut",distcut,"rep_cut",max_replicate_th,date,sep="_"))
    output_dir <- file.path(allout_dir,paste0(subject_name,"_",file_filter,"_",date,"_R"))

    ## run for subject, collect output from this function
    subject_output_list <- subject(subject_name,pdir,output_dir,sample_exclusion_th,max_replicate_th,sel_normal_sample,all_normal_samples)

    if (sample_name_change) {  
        phylogenetic_tree_with_boostrap_sampleNameChange(subject_name,pdir,output_dir,sample_exclusion_th,max_replicate_th,sel_normal_sample,all_normal_samples,colorfile_header,KLM,allout_dir,bscut,method='UPGMA')
        phylogenetic_tree_with_boostrap_sampleNameChange(subject_name,pdir,output_dir,sample_exclusion_th,max_replicate_th,sel_normal_sample,all_normal_samples,colorfile_header,KLM,allout_dir,bscut,method='NJ')
    } else {
        #phylogenetic_tree_with_boostrap()
        #phylogenetic_tree_with_boostrap_UPGMA()
    }

}



