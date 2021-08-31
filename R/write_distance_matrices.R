##' write_distance_matrices
##' @export
write_distance_matrices <- function(output_dir,fn_suffix,sample_names,marker_names,raw_dist,repre_replicates) {
    ### Generate distance matrix per marker for the selected replicates

    n <- length(sample_names)
    for (marker in marker_names) {
        results <- matrix(NA, nrow = n, ncol = n)

        for (s1 in 1:n) {
            sample1 <- sample_names[s1]
            rep_sample1 <- repre_replicates[[marker]][[sample1]]

            for (s2 in 1:n) {
                sample2 <- sample_names[s2]
                rep_sample2 <- repre_replicates[[marker]][[sample2]]

                if ((length(rep_sample1) > 0) & (length(rep_sample2) > 0)) {
                    results[s1,s2] <- raw_dist[[marker]][rep_sample1,rep_sample2]
                }
            }
        }
        ma_dist_df <- as.data.frame(results)
        colnames(ma_dist_df) <- sample_names
        rownames(ma_dist_df) <- sample_names
        write.table(ma_dist_df,file=file.path(output_dir,paste0(marker,"_",fn_suffix)),sep="\t",quote=FALSE,col.names=NA)
    }
}
