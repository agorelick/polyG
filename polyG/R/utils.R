
##' theme_std
##' @export
theme_std <- function(base_size = 11, base_line_size = base_size/22, base_rect_size = base_size/22) {
    base_fam <- 'ArialMT'
    base_face <- 'plain'
    require(ggplot2)
    theme_classic(base_size = base_size, base_family = base_fam,
                  base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace%
    theme(
          line = element_line(colour = "black", size = base_line_size, linetype = 1, lineend = "round"),
          text = element_text(family = base_fam, #face = base_face,
                    colour = "black", size = base_size, lineheight = 0.9,
                    hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug=F),
          axis.text = element_text(colour = "black", family=base_fam, #face=base_face, 
                                   size=rel(0.8)),
          axis.ticks = element_line(colour = "black", size=rel(1)),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", size = rel(1)),
          legend.key = element_blank(),
          strip.background = element_blank())
}



##' get_output_filepath
##' @export
get_output_filepath <- function(subject_name,title,suffix,output_dir,max_replicate_th) {
    fn <- paste0(subject_name,"_",title,"_",max_replicate_th,suffix)
    fn <- file.path(output_dir,fn)
    return(fn)             
}



##' read_distance_matrix
##' @export
read_distance_matrix <- function(file,return.as.matrix=T) {
    ## read txt file with a saved distance matrix (i.e. table with named rows and cols and numeric distances as cells)
    distance_matrix <- fread(file)
    rows <- distance_matrix[[1]] 
    distance_matrix <- distance_matrix[,(2:ncol(distance_matrix)),with=F]
    m <- as.matrix(distance_matrix)
    rownames(m) <- rows
    if(return.as.matrix==F) {
        as.dist(m,diag=T)
    } else {
        m
    }
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

