##' plot_overlap_scores.R
##' @export
plot_overlap_scores <- function(overlap_score,all_sample_names,marker_names,ps_filepath,psplot_filepath) {
    ### Plot overlap scores for all markers per sample
    require(ggplot2)

    Sample <- character()
    Overlap_Score <- numeric()
    nmarkers <- length(marker_names) 
    nsamples <- length(all_sample_names)

    for (sample in all_sample_names) {
        for (marker in marker_names) {
            Overlap_Score <- c(Overlap_Score,overlap_score[[marker]][sample])
        }
        Sample <- c(Sample,rep(sample,nmarkers))
    }
    Marker <- c(rep(marker_names,nsamples))
    Sample <- factor(Sample, levels=all_sample_names)
    df <- data.frame(Sample = Sample, Marker = Marker, Overlap_Score = Overlap_Score)
    write.table(df,file=ps_filepath,sep="\t",quote=FALSE,row.names=FALSE)

    p <- ggplot(df, aes(x=Sample, y=Overlap_Score))+
        geom_boxplot(aes(fill=Sample),outlier.colour = NA) + guides(fill='none') + 
        geom_jitter()+
        theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1))+
        theme(axis.text.y = element_text(size=12))+
        theme(axis.title.x = element_text(size=15,margin = margin(t = 20, r = 0, b = 0, l = 0)))+
        theme(axis.title.y = element_text(size=15,margin = margin(t = 0, r = 20, b = 0, l = 0))) #+
    ggsave(psplot_filepath, width=15, height=10, units='cm', dev=cairo_pdf)
}
