##' create_distance_matrix
##' @export
create_distance_matrix <- function(sample_names,final_marker_names,repre_replicates,raw_dist) {
    ### Create JSM distance matrix for representative replicates
    ### Only markers with no missing samples are used

    n <- length(sample_names)
    dm <- matrix(0, nrow = n, ncol = n)
    nmarkers <- matrix(0, nrow = n, ncol = n)
    dm_df <- as.data.frame(dm)
    colnames(dm_df) <- sample_names
    rownames(dm_df) <- sample_names

    for (marker in final_marker_names[["used_markers"]]) {
        for (s1 in 1:n) {
            sample1 <- sample_names[s1]
            rep_sample1 <- repre_replicates[[marker]][[sample1]]

            for (s2 in 1:n) {
                sample2 <- sample_names[s2]
                rep_sample2 <- repre_replicates[[marker]][[sample2]]

                if ((s2 != s1) & (length(rep_sample1) > 0) & (length(rep_sample2) > 0)) {
                    dm_df[s1,s2] <- dm_df[s1,s2] + raw_dist[[marker]][rep_sample1,rep_sample2]
                    nmarkers[s1,s2] <- nmarkers[s1,s2] + 1
                }
            }
        }
    }

    # normalize the pairwise distance by the number of compared markers
    for (s1 in 1:n) {
        for (s2 in 1:n) {
            if (nmarkers[s1,s2] > 0) {
                dm_df[s1,s2] <- dm_df[s1,s2]/nmarkers[s1,s2]
            }
        }
    }
    return(dm_df)
}
