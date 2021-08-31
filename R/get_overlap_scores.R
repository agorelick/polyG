##' get_overlap_scores
##' @export
get_overlap_scores <- function(prob_distr,all_sample_names,sel_normal_sample) {
    ### Calculate the difference in the polyG probability distribution between the
    ### representative replicate of a sample and that of the normal (reference) sample.
    ### Get the sum of the positive components of the difference.
    ### norm_df_repre_repli[[marker]]: normalized probability distributions of the  
    ###                                input polyg data; representative replicates only
    ### return: overlap score defined as "sum of the positive components of the difference
    ###         between the probability distribution of a sample and that of the normal"

    oscore <- numeric()
    for (sample in all_sample_names) {
        if (sample %in% colnames(prob_distr)) {
            sampidx <- which(colnames(prob_distr)==sample)
            nidx <- which(colnames(prob_distr)==sel_normal_sample)                
            prob_distr_diff <- prob_distr[,sampidx]-prob_distr[,nidx]
            oscore[sample] <- sum(prob_distr_diff[which(prob_distr_diff > 0)])
        } else {
            oscore[sample] <- NA
        }
    }
    return(oscore)
}
