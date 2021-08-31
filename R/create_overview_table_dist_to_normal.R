##' create_overview_table_dist_to_normal
##' @export
create_overview_table_dist_to_normal <- function(filepath,repre_replicates,raw_dist,all_sample_names,marker_names,sel_normal_sample) {
    ### Create file with overview about selected replicates and their distance to the normal
    Sample <- all_sample_names
    nsamples <- length(Sample)

    df <- data.frame(Sample, stringsAsFactors = F)
    rownames(df) <- all_sample_names

    df[,2] <- rep(0,nsamples)
    names(df)[2] <- paste0("JSDistance_",sel_normal_sample)

    i <- 3
    for (marker in marker_names) {
        for (sample in all_sample_names) {
            if ((length(repre_replicates[[marker]][[sample]]) != 0) & (length(repre_replicates[[marker]][[sel_normal_sample]]) != 0)) {
                df[sample,i] <- repre_replicates[[marker]][[sample]][1]  # name of selected replicate of marker
                df[sample,i+1] <- raw_dist[[marker]][repre_replicates[[marker]][[sample]][1],repre_replicates[[marker]][[sel_normal_sample]][1]]
                df[sample,2] <- df[sample,2] + df[sample,i+1]
            } else {
                df[sample,i] <- NA
                df[sample,i+1] <- NA
            }
        }
        names(df)[i] <- paste0(marker,"_selected")
        names(df)[i+1] <- paste0(marker,"_JSD")

        i <- i+2
    }
    write.table(df,file=filepath,sep="\t",row.names=F,col.names=T,quote=F)
}
