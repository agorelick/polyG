##' distance_distribution_heatmap
##' 
##' Perform average/UPGMA linkage on the JSM distance matrix and create distrance distribution heatmap. use final_marker_names and sample_names with possible exclusion of certain samples and markers
##' 
##' @export
distance_distribution_heatmap <- function(filepath_old,filepath_new,subject_name,repre_replicates,track_lengths,final_marker_names,markers,max_replicate_th,sel_normal_sample,all_normal_samples,sample_names,new_rownames,rob_th,add_info_dir,raw_dist) {
  nm <- length(final_marker_names[[markers]])
  ns <- length(sample_names)
  # no heatmap file will be generated if no usable markers are available
  if ((nm < 2) | (ns < 2) ) {
    print(paste0("Heatmap for ",markers," not generated. Need at least two usable markers and two usable samples."))
    
  } else {
    dis_distr <- matrix(NA,nrow=ns,ncol=nm)
    
    notecol <- rep("black",ns*nm)

    midx <- 0
    for (marker in final_marker_names[[markers]]) {
      
      # record ambigous indels for each marker
      ambiguous_indels_file <- file.path(add_info_dir,paste0(marker,"_",markers,"_ambiguous_indels.txt"))
      if (file.exists(ambiguous_indels_file)) {
        file.remove(ambiguous_indels_file)
      }
      
      midx <- midx + 1
      normal_sample <- character()
      if (length(repre_replicates[[marker]][[sel_normal_sample]]) > 0) {
        normal_sample <- sel_normal_sample
      } else {
        for (sample in all_normal_samples) {
          if (length(repre_replicates[[marker]][[sample]]) > 0) {
            normal_sample <- sample
            break
          }
        }
      }
      if (length(normal_sample) == 0) {
        dis_distr[,midx] <- rep(0,ns)
        print(paste0("Normal sample for marker ",marker," not available. No deletion/insertion classification!"))
      } else {  
        rep_normal <- repre_replicates[[marker]][[normal_sample]]
        lmean_normal <- track_lengths[[marker]][[rep_normal]]$lmean
        lmean3_normal <- track_lengths[[marker]][[rep_normal]]$lmean3
        lpeak_normal <- track_lengths[[marker]][[rep_normal]]$lpeak
        lpeak2_normal <- track_lengths[[marker]][[rep_normal]]$lpeak2
        lmedian_normal <- track_lengths[[marker]][[rep_normal]]$lmedian
        
        write(paste0(marker," in normal sample ",normal_sample," mean ",lmean_normal," mean3 ",lmean3_normal," peak ",lpeak_normal," peak2 ",lpeak2_normal," median ",lmedian_normal),file=ambiguous_indels_file,append=TRUE)
        
        sidx <- 0
        for (sample in sample_names) {
          sidx <- sidx + 1
          rep_sample <- repre_replicates[[marker]][[sample]]
          if (length(rep_sample) > 0) {
            if (sample == normal_sample) {
              dis_distr[sidx,midx] <- 0
              next
            }
            
            lmean <- track_lengths[[marker]][[rep_sample]]$lmean
            lmean3 <- track_lengths[[marker]][[rep_sample]]$lmean3
            lpeak <- track_lengths[[marker]][[rep_sample]]$lpeak
            lpeak2 <- track_lengths[[marker]][[rep_sample]]$lpeak2
            lmedian <- track_lengths[[marker]][[rep_sample]]$lmedian
            
            # likely deletion
            if (((lmean_normal > lmean) & (lmean3_normal > lmean3)) | 
                ((lmean_normal + lmean3_normal - rob_th) > (lmean + lmean3))) {
              
              dis_distr[sidx,midx] <- -raw_dist[[marker]][rep_normal,rep_sample]
              
            } else if (((lmean_normal < lmean) & (lmean3_normal < lmean3)) | 
                        ((lmean_normal + lmean3_normal + rob_th) < (lmean + lmean3))) {
              # likely insertion 
              
              dis_distr[sidx,midx] <- raw_dist[[marker]][rep_normal,rep_sample]
              
       #     } else if (lmean3 < lmean3_normal) {
       #     I think I would not agree with Hannes's choice of lmean3 here; lmean seems better.       
            } else if (lmean < lmean_normal) {
              # possibly deletion
              dis_distr[sidx,midx] <- -raw_dist[[marker]][rep_normal,rep_sample]
              write(paste0(markers,": ",marker," in sample ",sample," ambiguous deletion: mean ",lmean,
                           " mean3 ",lmean3," peak ",lpeak," peak2 ",lpeak2," median ",lmedian),
                    file=ambiguous_indels_file,append=TRUE)
            } else if (lmean > lmean_normal) {
              # possibly insertion
              dis_distr[sidx,midx] <- raw_dist[[marker]][rep_normal,rep_sample]
              write(paste0(markers,": ",marker," in sample ",sample," ambiguous insertion: mean ",lmean,
                           " mean3 ",lmean3," peak ",lpeak," peak2 ",lpeak2," median ",lmedian),
                    file=ambiguous_indels_file,append=TRUE)
            } else {
              # unable to classify sample as insertion or deletion
              dis_distr[sidx,midx] <- -raw_dist[[marker]][rep_normal,rep_sample]
              write(paste0(markers,": ",marker," in sample ",sample," Could not classify insertion or deletion: mean ",lmean,
                           " mean3 ",lmean3," peak ",lpeak," peak2 ",lpeak2," median ",lmedian),
                    file=ambiguous_indels_file,append=TRUE)
              
            }
          }
        }
      }
    }
    
    #browser()
    short_marker_names <- sapply(final_marker_names[[markers]],function(x){unlist(strsplit(x,paste0("_",subject_name)))[1]})
    names(short_marker_names) <- NULL
    
    dis_distr_df <- as.data.frame(dis_distr)
    
    ## heatmap with original sample names
    colnames(dis_distr_df) <- short_marker_names
    rownames(dis_distr_df) <- sample_names
    #rownames(dis_distr_df) <- new_rownames
    
    dis_distr <- data.matrix(dis_distr_df)
    dis_distr <- round(dis_distr,2)
    color <- colorRampPalette(c("seagreen", "white", "purple3"))
    
    pdf(filepath_old,width=10,height=7,pointsize = 9)
    par(mar=c(10,4,4,2))
    ## for RDS
    #heatmap.2(dis_distr,trace="none",na.rm=FALSE,cellnote=dis_distr,
    #          notecex=1/log10(0.61*nm),notecol="black",na.color="yellow",cexRow = max(1.8,0.3+1/log10(ns)), 
    #          notecex=1/log10(0.61*nm),notecol=notecol,na.color="yellow",cexRow = max(1.8,0.3+1/log10(ns)), 
    #          cexCol = max(2,0.3+1/log10(nm)),density.info="none",col=color,
    #          lhei=c(1,6),lwid=c(1,6),sepwidth=c(0.005,0.005),sepcolor="white",
    #          colsep=1:nm,rowsep=1:ns,margins = c(9,7))
    # for Emma's data
    heatmap.2(dis_distr,trace="none",na.rm=FALSE,cellnote=dis_distr,
              notecex=0.01+1/(1.6*log10(nm)),notecol=notecol,na.color="yellow",cexRow = min(1,0.3+1/log10(ns)), 
              cexCol = min(1.2,0.3+1/log10(nm)),density.info="none",col=color,
              lhei=c(1,6),lwid=c(1,4),sepwidth=c(0.005,0.005),sepcolor="white",
              colsep=1:nm,rowsep=1:ns,margins = c(5,6.5))
    ## for data with few samples as in UW and Kosima  
    #pdf(filepath,pointsize = 9,width=6,height=4)
    #par(mar=c(10,4,4,2))
    #heatmap.2(dis_distr,trace="none",na.rm=FALSE,cellnote=dis_distr,
    #          notecex=0.2+1/(1.2*log10(nm)),notecol=notecol,na.color="yellow",cexRow = min(1.1,0.6+1/log10(ns)), 
    #          cexCol = min(1,0.3+1/log10(nm)),density.info="none",col=color,
    #          lhei=c(1,3),lwid=c(1,4),sepwidth=c(0.005,0.005),sepcolor="white",
    #          colsep=1:nm,rowsep=1:ns,margins = c(8,8),srtCol=45)
    
    dev.off()
    
    ## heatmap with new sample names
    colnames(dis_distr_df) <- short_marker_names
    #rownames(dis_distr_df) <- sample_names
    rownames(dis_distr_df) <- new_rownames
    
    dis_distr <- data.matrix(dis_distr_df)
    dis_distr <- round(dis_distr,2)
    color <- colorRampPalette(c("seagreen", "white", "purple3"))
    
    pdf(filepath_new,width=10,height=7,pointsize = 9)
    par(mar=c(10,4,4,2))
    ## for RDS
    #heatmap.2(dis_distr,trace="none",na.rm=FALSE,cellnote=dis_distr,
    #          #notecex=1/log10(0.61*nm),notecol="black",na.color="yellow",cexRow = max(2,0.3+1/log10(ns)), 
    #          notecex=1/log10(0.61*nm),notecol=notecol,na.color="yellow",cexRow = max(2,0.3+1/log10(ns)), 
    #          cexCol = max(2,0.3+1/log10(nm)),density.info="none",col=color,
    #          lhei=c(1,6),lwid=c(1,6),sepwidth=c(0.005,0.005),sepcolor="white",
    #          colsep=1:nm,rowsep=1:ns,margins = c(9,6.5))
    
    # for Emma's data
    heatmap.2(dis_distr,trace="none",na.rm=FALSE,cellnote=dis_distr,
              notecex=0.01+1/(1.6*log10(nm)),notecol=notecol,na.color="yellow",cexRow = min(1,0.3+1/log10(ns)), 
              cexCol = min(1.2,0.3+1/log10(nm)),density.info="none",col=color,
              lhei=c(1,6),lwid=c(1,4),sepwidth=c(0.005,0.005),sepcolor="white",
              colsep=1:nm,rowsep=1:ns,margins = c(5,6.5))
    
    dev.off()
    
  }
}
