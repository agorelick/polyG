##' distance_distribution_plot
##' @export
distance_distribution_plot <- function(dis_distr_fp,repre_replicates,track_lengths,raw_dist,all_sample_names,marker_names,sel_normal_sample,max_replicate_th) {
    ### Visualize the distances of all samples of a marker to the selected normal sample to get an idea of the evolutionary information provided 

    require("colorspace")
    m <- length(repre_replicates)
    n <- length(all_sample_names)
    dis_distr <- matrix(NA, nrow = n-1, ncol = m)

    midx <- 1
    for (marker in marker_names) {
        sidx <- 1
        for (sample in all_sample_names) {
            if (sample != sel_normal_sample) {
                sample_rep <- repre_replicates[[marker]][[sample]]
                normal_rep <- repre_replicates[[marker]][[sel_normal_sample]]
                if ((length(sample_rep) > 0) & (length(normal_rep) > 0)) {
                    # likely insertion
                    if ((track_lengths[[marker]][[normal_rep]]$lmean < track_lengths[[marker]][[sample_rep]]$lmean) 
                        & (track_lengths[[marker]][[normal_rep]]$lmean3 < track_lengths[[marker]][[sample_rep]]$lmean3)) {
                        dis_distr[sidx,midx] <- raw_dist[[marker]][sample_rep,normal_rep]
                    } 
                    # likely deletion
                    if ((track_lengths[[marker]][[normal_rep]]$lmean > track_lengths[[marker]][[sample_rep]]$lmean)
                        & (track_lengths[[marker]][[normal_rep]]$lmean3 > track_lengths[[marker]][[sample_rep]]$lmean3)) { 
                        dis_distr[sidx,midx] <- -raw_dist[[marker]][sample_rep,normal_rep]  
                    }
                }
                sidx <- sidx + 1
            }
        }
        midx <- midx + 1
    }
    dis_distr_df <- as.data.frame(dis_distr)
    colnames(dis_distr_df) <- marker_names

    color <- rainbow_hcl(n=length(marker_names), c=80, l=60, 0, 300)
    pdf(dis_distr_fp)
    boxplot(dis_distr_df,names=marker_names,col=color,xlab="Marker",ylab=paste0("Jensen-Shannon distance to normal sample ",sel_normal_sample),cex.lab=1.4)
    abline(h=max_replicate_th,col="red",lty=3,lwd=3)
    abline(h=-max_replicate_th,col="red",lty=3,lwd=3)

    dev.off()
}
