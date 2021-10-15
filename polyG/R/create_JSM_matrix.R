##' create_JSM_matrix
##' @export
create_JSM_matrix <- function(norm_df, filepath) {

    ### Calculate Jensen-Shannon distance for each pair of samples of the given marker
    ### and creates a matrix with the pairwise distances
    ### norm_df: normalized probability distributions of the input polyg data
    ### return: Jensen-Shannon distance matrix

    require(phyloseq)

    jsm_df <- as.matrix(sqrt((distance(otu_table(norm_df,taxa_are_rows = T),"jsd"))/log(2)))
    colnames(jsm_df) <- colnames(norm_df)
    rownames(jsm_df) <- colnames(norm_df)
    write.table(jsm_df,file=filepath,sep="\t",quote=FALSE,col.names=NA)
    return(jsm_df)
}

