##' angular_distance
##'
##' Generate angular distance matrices, heatmaps, and trees. Optionally generate bootstrapped data.
##'
##' @export
angular_distance <- function(input_dir, allout_dir, subject_name, sel_normal_sample, all_normal_samples, new_sample_names, bootstrap, bscut) { 

    # hardcoded parameters
    markergroup <- "usedmarkers"
    replicate <- "repreReplicate"
    power <- 1

    # output from subject()
    output_dir <- file.path(allout_dir,paste0(subject_name,"_R"))
    
    # output directories for angular distance
    angular_dist_dir <- file.path(allout_dir,"results_angular_distance_representativeReplicates")
    dir.create(angular_dist_dir,recursive = TRUE, showWarnings = FALSE)
    ad_matrix_w_root_dir <- file.path(angular_dist_dir,"angular_dist_matrix_w_root_usedmarkers")
    dir.create(ad_matrix_w_root_dir,recursive = TRUE, showWarnings = FALSE) 
    ad_heatmap_dir <- file.path(angular_dist_dir,paste0("angular_dist_heatmap_",markergroup))
    dir.create(ad_heatmap_dir,recursive = TRUE, showWarnings = FALSE)  
    ad_tree_w_root_dir <- file.path(angular_dist_dir,"angular_dist_trees_w_root_usedmarkers")
    dir.create(ad_tree_w_root_dir,recursive = TRUE, showWarnings = FALSE)

    #* only use "used markers" and non-excluded samples for trees
    # get the list of used markers
    addinfodir <- file.path(output_dir,"additional_info")
    mn_info_file <- grep('final_marker_names.txt',dir(addinfodir,full.names=T),value=T)
    mn_info <- read.table(mn_info_file,sep=" ",stringsAsFactors = FALSE,check.names = FALSE)
    um_row <- grep("used_markers",mn_info[,1])
    um_long <- mn_info[um_row,2]
    usedmarkers <- sapply(um_long,function(x) unlist(strsplit(x,"_"))[1])
    nusedmarkers <- length(usedmarkers)

    #* get normalized_data 
    normdatadir <- file.path(output_dir,"normalized_data")
    normdf_files <- dir(normdatadir,full.names=T)

    # all input samples
    samplenames_file <- file.path(input_dir,"sample_names",paste0(subject_name,"_SampleNames.txt"))
    allsamps <- (read.table(samplenames_file,sep="\t",stringsAsFactors=FALSE))$V1

    # color file (sample colors for trees including multiple annotations)
    colorfile <- read.table(file.path(input_dir,"color_files",paste0(subject_name,"_colorfile.txt")),sep="\t",header=T,stringsAsFactors=FALSE,fill=TRUE)

    meanlen <- list()    
    meanlen_minus_normal <- list()  
    exclude_samps <- character() 
    for (marker_file in normdf_files) {
        marker <- tail(strsplit(marker_file,'[/]')[[1]],1)
        marker <- strsplit(marker,"_")[[1]][1]
        if (!(marker %in% usedmarkers)) next

        normdf_repre_repli <- read.table(marker_file,sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
        meanlen[[marker]] <- rep(0,ncol(normdf_repre_repli))
        for (row in 1:nrow(normdf_repre_repli)) {
            row_frac <- as.numeric(normdf_repre_repli[row,])
            meanlen[[marker]] <- meanlen[[marker]] + row*row_frac
        }

        samp <- colnames(normdf_repre_repli)
        col_norm <- which(samp %in% all_normal_samples)
        meanlen_normals <- mean(meanlen[[marker]][col_norm])
        meanlen_minus_normal[[marker]] <- meanlen[[marker]] - meanlen_normals
        names(meanlen_minus_normal[[marker]]) <- samp
        sel_excl <- !(allsamps %in% samp)
        exclude_samps <- c(exclude_samps, allsamps[sel_excl])
    }
    exclude_samps <- unique(exclude_samps)  
    excl <- c(all_normal_samples,exclude_samps)
    usedsamps <- allsamps[-which(allsamps %in% excl)]
    ns <- length(usedsamps)

    #* get Z
    meanlen_diff_normal <- meanlen_minus_normal
    markerset <- usedmarkers
    Z <- Zij(usedsamps,markerset,meanlen_diff_normal)

    #* get heatmap - use normal as reference
    usedsamps_nor <- c(usedsamps,sel_normal_sample)
    new_rownames <- usedsamps_nor
    if (sample_name_change) {
        new_rownames <- new_sample_names[pmatch(usedsamps_nor,new_sample_names[,1]),2]
    } 

    if (ns < 1) {
        cat(paste0(subject_name," ",markergroup," ",replicate," has fewer than 2 samples available. No angular distance heatmap will be produced. \n"))
        write(paste0(subject_name," ",markergroup," ",replicate," has fewer than 2 samples available. No angular distance heatmap will be produced. \n"),file=file.path(ad_heatmap_dir,paste0(subject_name,"_",markergroup,"_",replicate,"_error.txt")))
    } else {
        adh_table_fp_old <- file.path(ad_heatmap_dir,paste0(subject_name,"_angular_dist_heatmap_table_",markergroup,"_oldnames.txt"))
        adh_table_fp_new <- file.path(ad_heatmap_dir,paste0(subject_name,"_angular_dist_heatmap_table_",markergroup,"_newnames.txt"))
        adh_fp_old <- file.path(ad_heatmap_dir,paste0(subject_name,"_angular_dist_heatmap_",markergroup,"_oldnames.pdf"))
        adh_fp_new <- file.path(ad_heatmap_dir,paste0(subject_name,"_angular_dist_heatmap_",markergroup,"_newnames.pdf"))
        angular_distance_heatmap(adh_fp_old,adh_fp_new,adh_table_fp_old,adh_table_fp_new,subject_name,Z,markerset,usedsamps,sel_normal_sample,new_rownames,sample_name_change, markergroup, replicate)
    }

    #* get angular distance matrix
    if (ns < 2) {
        cat(paste0(subject_name," ",markergroup," ",replicate," has fewer than 2 samples available. No angular distance matrix will be produced. \n"))
        write(paste0(subject_name," ",markergroup," ",replicate," has fewer than 2 samples available. No angular distance matrix will be produced. \n"),file=file.path(ad_matrix_w_root_dir,paste0(subject_name,"_",markergroup,"_",replicate,"_error.txt")))
    } else {
        angular_dist_w_root <- ang_dist_matrix(Z,usedsamps,ns,markerset,sel_normal_sample,power)
        write.table(angular_dist_w_root,file=file.path(ad_matrix_w_root_dir,paste0(subject_name,"_angular_dist_matrix_w_root_",markergroup,"_",replicate,"_oldnames.txt")),sep="\t",quote=FALSE,col.names=NA)

        if (sample_name_change) {
            angular_dist_w_root_newnames <- angular_dist_w_root
            colnames(angular_dist_w_root_newnames) <- new_sample_names[pmatch(colnames(angular_dist_w_root),new_sample_names[,1]),2]
            rownames(angular_dist_w_root_newnames) <- colnames(angular_dist_w_root_newnames)
            write.table(angular_dist_w_root_newnames,file=file.path(ad_matrix_w_root_dir,paste0(subject_name,"_angular_dist_matrix_w_root_",markergroup,"_",replicate,"_newnames.txt")),sep="\t",quote=FALSE,col.names=NA)
        }
    }

    ##* Build phylogenetic trees with the angular distance 
    if (ns < 2) {
        cat(paste0(subject_name," ",markergroup," ",replicate," has fewer than 2 samples available. No angular distance tree will be produced. \n"))
        write(paste0(subject_name," ",markergroup," ",replicate," has fewer than 2 samples available. No angular distance tree will be produced. \n"),file=file.path(ad_tree_w_root_dir,paste0(subject_name,"_",markergroup,"_",replicate,"_error.txt")))
    } else {
        generate_angular_distance_trees(angular_dist_w_root, sel_normal_sample, colorfile, sample_name_change, ad_tree_w_root_dir, markerset, usedsamps, ns, power, new_sample_names, meanlen_diff_normal, bootstrap, bscut, subject_name, markergroup, replicate)  
    }
}


