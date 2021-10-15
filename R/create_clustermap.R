##' create_clustermap
##' @export
create_clustermap <- function(output_dir,subject_name,dm_df,max_replicate_th,sample_names,anno,sample_name_change) {
### Perform hierarchical clustering on the normalized JSM distance matrix and create the clustermap

require(gplots)
require(RColorBrewer)
require(colorspace)

  dm_filepath <- get_output_filepath(subject_name,"distancematrix",".tsv",output_dir,max_replicate_th)
  write_distance_matrix(dm_filepath,dm_df)
  
  ns <- length(sample_names)
  
  dm <- data.matrix(dm_df)
  dm_round <- round(dm,2)
  
  color <- colorRampPalette(c("white", "rosybrown2", "violetred4"))
  
  if (!sample_name_change) {
    cm_filepath = get_output_filepath(subject_name,"clustermap",".pdf",output_dir,max_replicate_th)
    pdf(cm_filepath,pointsize = 7)
    
    heatmap.2(dm,Rowv=TRUE,symm=TRUE,dendrogram="row",trace="none",na.rm=FALSE,cellnote=dm_round,notecex=0.8,notecol="black",na.color="yellow",cexRow = 0.3+1/log10(ns), cexCol = 0.3+1/log10(ns),density.info="none",col=color,lhei=c(1,8),lwid=c(0.5,3.5),margins=c(8,8))
    dev.off()
  } else {
    # cluster map with old sample names
    cm_filepath_old = get_output_filepath(subject_name,"clustermap_oldnames",".pdf",output_dir,max_replicate_th)
    
    pdf(cm_filepath_old,pointsize = 7)
    heatmap.2(dm,Rowv=TRUE,symm=TRUE,dendrogram="row",trace="none",na.rm=FALSE,cellnote=dm_round,notecex=0.8,notecol="black",na.color="yellow",cexRow = 0.3+1/log10(ns), cexCol = 0.3+1/log10(ns),density.info="none",col=color,lhei=c(1,8),lwid=c(0.5,3.5),margins=c(8,8))
    dev.off()
    
    # cluster map with new sample names
    dm_new <- dm
    colnames(dm_new) <- anno[pmatch(colnames(dm),anno[,1]),2]
    rownames(dm_new) <- anno[pmatch(rownames(dm),anno[,1]),2]
    
    cm_filepath_new = get_output_filepath(subject_name,"clustermap_newnames",".pdf",output_dir,max_replicate_th)
    
    pdf(cm_filepath_new,pointsize = 7)
    heatmap.2(dm_new,Rowv=TRUE,symm=TRUE,dendrogram="row",trace="none",na.rm=FALSE,cellnote=dm_round,notecex=0.8,notecol="black",na.color="yellow",cexRow = 0.3+1/log10(ns), cexCol = 0.3+1/log10(ns),density.info="none",col=color,lhei=c(1,8),lwid=c(0.5,3.5),margins=c(8,8))
    dev.off()
  }
  
}
