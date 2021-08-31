##' get_output_filepath
##' @export
get_output_filepath <- function(subject_name,title,suffix,output_dir,max_replicate_th) {
    fn <- paste0(subject_name,"_",title,"_",max_replicate_th,suffix)
    fn <- file.path(output_dir,fn)
    return(fn)             
}

##' write_distance_matrix
##' @export
write_distance_matrix <- function(filepath,dm_df) {
    write.table(dm_df,file=filepath,sep="\t",quote=FALSE,col.names=NA)
}

##' create_dir_check_if_already_exists
##' @export
create_dir_check_if_already_exists <- function(dir,overwrite) {
    ## Check if the expected output directories already exist. If so, prompt user to specify overwrite argument and throw error.
    if(dir.exists(dir) & overwrite==F) {
        stop(paste(dir,'already exists! To re-run, delete this directory or add argument overwrite=TRUE.')) 
    } else if(dir.exists(dir) & overwrite==T){
        message(paste(dir,'already exists! proceeding anyway.'))
    } else if(!dir.exists(dir)) {
        message('Created directory: ',dir)
        dir.create(dir,recursive=T,showWarnings=F)
    }
}

##' error_if_file_not_exist
##' @export
error_if_file_not_exist <- function(file) {
    ## Check if the expected input files exist. If not, prompt user to fix and throw error.
    if(!file.exists(file)) stop(paste(file,'does not exist!'))
}

