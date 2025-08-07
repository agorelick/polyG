##' get_output_filepath
##'
##' Generate filename and path
##'
##' @export
get_output_filepath <- function(subject_name,title,suffix,output_dir,max_replicate_th) {
    fn <- paste0(subject_name,"_",title,"_",max_replicate_th,suffix)
    fn <- file.path(output_dir,fn)
    return(fn)             
}


##' write_distance_matrix
##'
##' Utilty to write distance matrices in standard format
##'
##' @export
write_distance_matrix <- function(filepath,dm_df) {
    write.table(dm_df,file=filepath,sep="\t",quote=FALSE,col.names=NA)
}


##' create_JSM_matrix
##'
##' Calculate Jensen-Shannon distance for each pair of samples of the given marker and creates a matrix with the pairwise distances norm_df: normalized probability distributions of the input polyg data. Returns: Jensen-Shannon distance matrix
##' 
##' @export
create_JSM_matrix <- function(norm_df) {
    jsm_df <- as.matrix(sqrt((distance(otu_table(norm_df,taxa_are_rows = TRUE),"jsd"))/log(2)))
    colnames(jsm_df) <- colnames(norm_df)
    rownames(jsm_df) <- colnames(norm_df)
    return(jsm_df)
}


##' create_distance_matrix
##' 
##' Create JSM distance matrix for representative replicate. Only markers with no missing samples are used. (Deprecated?)
##' 
##' @export
create_distance_matrix <- function(sample_names,final_marker_names,repre_replicates,raw_dist) {
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

##' select_representative_replicates
##' 
##' Select most representative replicates. First check if all of them agree, then select the pair which minimizes the provided distance and then select the replicate that maximizes the intensity
##' 
##' @export
select_representative_replicates <- function(marker,raw_dist,intensity,sample_replicates,repre_replicates,max_replicate_th,all_sample_names,outlist) {
    problematic_samples <- character()  # samples whose replicates are too dissimilar
    missing_samples <- character()  # samples missing from marker file for a given marker

    for (sample in all_sample_names) {
        if (length(sample_replicates[[marker]][[sample]]) == 0) {
            missing_samples <- append(missing_samples,sample)

        } else if (length(sample_replicates[[marker]][[sample]]) == 1) {
            repre_replicates[[marker]][[sample]] <- sample_replicates[[marker]][[sample]][1]

        } else if (length(sample_replicates[[marker]][[sample]]) == 2) {
            if (raw_dist[[marker]][sample_replicates[[marker]][[sample]][1],sample_replicates[[marker]][[sample]][2]] > max_replicate_th) {
                problematic_samples <- append(problematic_samples,sample)
            } else {
                int1 <- intensity[sample_replicates[[marker]][[sample]][1]]
                int2 <- intensity[sample_replicates[[marker]][[sample]][2]]
                if (int1 > int2) {
                    repre_replicates[[marker]][[sample]] <- sample_replicates[[marker]][[sample]][1]
                } else {
                    repre_replicates[[marker]][[sample]] <- sample_replicates[[marker]][[sample]][2]
                }
            }

        } else if (length(sample_replicates[[marker]][[sample]]) == 3) {
            repli <- sample_replicates[[marker]][[sample]]
            int <- intensity[repli]

            distance <- c(raw_dist[[marker]][repli[1],repli[2]],raw_dist[[marker]][repli[2],repli[3]],raw_dist[[marker]][repli[3],repli[1]])
            if (all(distance > max_replicate_th)) {
                problematic_samples <- append(problematic_samples,sample)
            } else {
                if (length(which(distance == min(distance))) > 1) {
                    repre_idx <- which(intensity[repli]==max(intensity[repli]))[1]
                    repre_replicates[[marker]][[sample]] <- repli[repre_idx]
                    print(paste0("For ",marker,", ",sample," has more than one pair of replicates with the minimal distance. Check sizing window."))
                } else {
                    min_dist_idx <- which(distance == min(distance))
                    if (min_dist_idx == 3) {
                        next_idx <- 1
                    } else {
                        next_idx <- min_dist_idx + 1
                    }
                    if (int[min_dist_idx] > int[next_idx]) {
                        repre_replicates[[marker]][[sample]] <- repli[min_dist_idx]
                    } else {
                        repre_replicates[[marker]][[sample]] <- repli[next_idx]
                    }
                }
            }

        } else if (length(sample_replicates[[marker]][[sample]]) == 4) {
            repli <- sample_replicates[[marker]][[sample]]
            int <- intensity[repli]

            distance <- c(raw_dist[[marker]][repli[1],repli[2]],raw_dist[[marker]][repli[2],repli[3]],raw_dist[[marker]][repli[3],repli[4]],raw_dist[[marker]][repli[4],repli[1]],raw_dist[[marker]][repli[1],repli[3]],raw_dist[[marker]][repli[2],repli[4]])
            if (all(distance > max_replicate_th)) {
                problematic_samples <- append(problematic_samples,sample)
            } else {
                min_dist_idx <- which(distance == min(distance))

                if (length(min_dist_idx) > 3) {
                    repre_idx <- which(intensity[repli]==max(intensity[repli]))[1]
                    repre_replicates[[marker]][[sample]] <- repli[repre_idx]
                    print(paste0("For ",marker,", ",sample," has more than three pairs of replicates with the minimal distance. Check sizing window."))

                } else if (length(min_dist_idx) == 3) {

                    if (setequal(min_dist_idx,c(1,2,5))) {
                        reps <- c(repli[1],repli[2],repli[3])
                        reps_idx <- which(intensity[reps]==max(intensity[reps]))[1]
                        repre_replicates[[marker]][[sample]] <- reps[reps_idx]
                        print(paste0("For ",marker,", ",sample," has three pairs of replicates with the minimal distance. Check sizing window."))

                    } else if (setequal(min_dist_idx,c(1,4,6))) {
                        reps <- c(repli[1],repli[2],repli[4])
                        reps_idx <- which(intensity[reps]==max(intensity[reps]))[1]
                        repre_replicates[[marker]][[sample]] <- reps[reps_idx]
                        print(paste0("For ",marker,", ",sample," has three pairs of replicates with the minimal distance. Check sizing window."))

                    } else if (setequal(min_dist_idx,c(2,3,6))) {
                        reps <- c(repli[2],repli[3],repli[4])
                        reps_idx <- which(intensity[reps]==max(intensity[reps]))[1]
                        repre_replicates[[marker]][[sample]] <- reps[reps_idx]
                        print(paste0("For ",marker,", ",sample," has three pairs of replicates with the minimal distance. Check sizing window."))

                    } else if (setequal(min_dist_idx,c(3,4,5))) {
                        reps <- c(repli[1],repli[3],repli[4])
                        reps_idx <- which(intensity[reps]==max(intensity[reps]))[1]
                        repre_replicates[[marker]][[sample]] <- reps[reps_idx]
                        print(paste0("For ",marker,", ",sample," has three pairs of replicates with the minimal distance. Check sizing window."))

                    } else {
                        repre_idx <- which(intensity[repli]==max(intensity[repli]))[1]
                        repre_replicates[[marker]][[sample]] <- repli[repre_idx]
                        print(paste0("For ",marker,", ",sample," has three pairs of replicates with the minimal distance. Check sizing window."))
                    }

                } else if (length(min_dist_idx) == 2) {

                    if (setequal(min_dist_idx,c(1,2)) | setequal(min_dist_idx,c(1,5)) | setequal(min_dist_idx,c(2,5))) {
                        reps <- c(repli[1],repli[2],repli[3])
                        reps_idx <- which(intensity[reps]==max(intensity[reps]))[1]
                        repre_replicates[[marker]][[sample]] <- reps[reps_idx]
                        print(paste0("For ",marker,", ",sample," has two pairs of replicates with the minimal distance. Check sizing window."))

                    } else if (setequal(min_dist_idx,c(1,4)) | setequal(min_dist_idx,c(1,6)) | setequal(min_dist_idx,c(4,6))) {
                        reps <- c(repli[1],repli[2],repli[4])
                        reps_idx <- which(intensity[reps]==max(intensity[reps]))[1]
                        repre_replicates[[marker]][[sample]] <- reps[reps_idx]
                        print(paste0("For ",marker,", ",sample," has two pairs of replicates with the minimal distance. Check sizing window."))

                    } else if (setequal(min_dist_idx,c(2,3)) | setequal(min_dist_idx,c(2,6)) | setequal(min_dist_idx,c(3,6))) {
                        reps <- c(repli[2],repli[3],repli[4])
                        reps_idx <- which(intensity[reps]==max(intensity[reps]))[1]
                        repre_replicates[[marker]][[sample]] <- reps[reps_idx]
                        print(paste0("For ",marker,", ",sample," has two pairs of replicates with the minimal distance. Check sizing window."))

                    } else if (setequal(min_dist_idx,c(3,4)) | setequal(min_dist_idx,c(3,5)) | setequal(min_dist_idx,c(4,5))) {
                        reps <- c(repli[1],repli[3],repli[4])
                        reps_idx <- which(intensity[reps]==max(intensity[reps]))[1]
                        repre_replicates[[marker]][[sample]] <- reps[reps_idx]
                        print(paste0("For ",marker,", ",sample," has two pairs of replicates with the minimal distance. Check sizing window."))

                    } else {
                        repre_idx <- which(intensity[repli]==max(intensity[repli]))[1]
                        repre_replicates[[marker]][[sample]] <- repli[repre_idx]
                        print(paste0("For ",marker,", ",sample," has two pairs of replicates with the minimal distance. Check sizing window."))
                    }

                } else {
                    if (min_dist_idx < 5) {
                        if (min_dist_idx == 4) {
                            next_idx <- 1
                        } else {
                            next_idx <- min_dist_idx + 1
                        }
                        if (int[min_dist_idx] > int[next_idx]) {
                            repre_replicates[[marker]][[sample]] <- repli[min_dist_idx]
                        } else {
                            repre_replicates[[marker]][[sample]] <- repli[next_idx]
                        }
                    } else if (min_dist_idx == 5) {
                        if (int[1] > int[3]) {
                            repre_replicates[[marker]][[sample]] <- repli[1]
                        } else {
                            repre_replicates[[marker]][[sample]] <- repli[3]
                        }
                    } else if (min_dist_idx == 6) {
                        if (int[2] > int[4]) {
                            repre_replicates[[marker]][[sample]] <- repli[2]
                        } else {
                            repre_replicates[[marker]][[sample]] <- repli[4]
                        }
                    }          
                }
            }
        }
    }

    nmiss <- length(missing_samples) + length(problematic_samples)
    all_missing_samples <- list()
    all_missing_samples[[marker]] <- union(missing_samples,problematic_samples)
    outlist[[marker]] <- list("all_missing_samples" = all_missing_samples[[marker]],"samples_not_in_marker_files" = missing_samples,"samples_wo_repre_repli" = problematic_samples,"repre_replicates" = repre_replicates[[marker]])
    return(outlist[[marker]])
}


##' Zij
##' 
##' function to calculate Zij and angular distance matrix
##' @export
Zij <- function(usedsamps,markerset,meanlen_diff_normal) {
    z <- list()
    denom <- numeric()
    for (t in usedsamps) {
        denom[t] <- 0
        for (m in markerset) {
            denom[t] <- denom[t] + (meanlen_diff_normal[[m]][t])^2
        }
        denom[t] <- sqrt(denom[t])
        for (m in markerset) {
            z[[m]][t] <- meanlen_diff_normal[[m]][t]/denom[t]
        }
    }
    z
}


##' ang_dist_matrix
##' 
##' generate angular distance matrix
##'
##' @export
ang_dist_matrix <- function(Z,usedsamps,ns,markerset,sel_normal_sample,power) {
    # for each pair of samples in Z, calculate angular distance (exclude excluded samples)
    # angular distance matrix
    angular_dist <- matrix(0,nrow=ns,ncol=ns)
    for (s1 in 1:(ns-1)) {
        for (s2 in (s1+1):ns) {
            for (m in markerset) {
                angular_dist[s1,s2] <- angular_dist[s1,s2] + Z[[m]][usedsamps[s1]]*Z[[m]][usedsamps[s2]]
            }
            angular_dist[s1,s2] <- (acos(pmin(pmax(angular_dist[s1,s2],-1.0),1.0)))^power
            angular_dist[s2,s1] <- angular_dist[s1,s2]
        }
    }
    colnames(angular_dist) <- usedsamps
    rownames(angular_dist) <- usedsamps

    #* add a column and row for the normal
    angular_dist_w_root <- rbind(angular_dist,rep((pi/3)^power,ns))
    angular_dist_w_root <- cbind(angular_dist_w_root,c(rep((pi/3)^power,ns),0))
    colnames(angular_dist_w_root) <- c(usedsamps,sel_normal_sample)
    rownames(angular_dist_w_root) <- c(usedsamps,sel_normal_sample)
    angular_dist_w_root
}


##' angular_distance_heatmap
##' 
##' Create heatmap for angular distrance w.r.t the normal 
##' 
##' @export
angular_distance_heatmap <- function(filepath_old,filepath_new,table_filepath_old,table_filepath_new,subject_name,Z,markerset,usedsamps,sel_normal_sample,new_rownames,sample_name_change,markergroup, replicate) {

    nm <- length(markerset)
    ns <- length(usedsamps)
    old_rownames <- c(usedsamps,sel_normal_sample)

    # no heatmap file will be generated if not enough usable markers or samples are available
    if ((nm < 2) | (ns < 1) ) {
        print(paste0("Angular distance heatmap for ",markergroup," ",replicate," not generated. Need at least two usable markers and two usable samples."))

    } else {
        dis_distr <- matrix(NA,nrow=(ns+1),ncol=nm)

        notecol <- rep("black",(ns+1)*nm)

        midx <- 0
        # For each marker, get the Z[[marker]][t] for each tumor sample
        for (marker in markerset) {
            midx <- midx + 1

            sidx <- 0
            for (t in old_rownames) {
                sidx <- sidx + 1
                if (t==sel_normal_sample) {
                    dis_distr[sidx,midx] <- 0
                } else {
                    dis_distr[sidx,midx] <- Z[[marker]][t] 
                }
            }
        }
    }

    dis_distr_df <- as.data.frame(dis_distr)

    ## heatmap with original sample names
    colnames(dis_distr_df) <- markerset
    rownames(dis_distr_df) <- old_rownames

    dis_distr <- data.matrix(dis_distr_df)
    write.table(dis_distr,file=table_filepath_old,sep="\t",quote=FALSE,col.names=NA)
    dis_distr <- round(dis_distr,2)
    color <- colorRampPalette(c("seagreen", "white", "purple3"))

    pdf(filepath_old,width=10,height=7,pointsize = 9)
    par(mar=c(10,4,4,2))

    heatmap.2(dis_distr,trace="none",na.rm=FALSE,cellnote=dis_distr,
              notecex=0.1+1/(1.6*log10(nm)),notecol=notecol,na.color="yellow",cexRow = min(1,0.3+1/log10(ns)), 
              cexCol = min(1.2,0.3+1/log10(nm)),density.info="none",col=color,
              lhei=c(1,6),lwid=c(1,4),sepwidth=c(0.005,0.005),sepcolor="white",
              colsep=1:nm,rowsep=1:ns,margins = c(5,6.5))
    dev.off()

    if (sample_name_change) {
        ## heatmap with new sample names
        colnames(dis_distr_df) <- markerset
        rownames(dis_distr_df) <- new_rownames

        dis_distr <- data.matrix(dis_distr_df)
        write.table(dis_distr,file=table_filepath_new,sep="\t",quote=FALSE,col.names=NA)
        dis_distr <- round(dis_distr,2)
        color <- colorRampPalette(c("seagreen", "white", "purple3"))

        pdf(filepath_new,width=10,height=7,pointsize = 9)
        par(mar=c(10,4,4,2))

        heatmap.2(dis_distr,trace="none",na.rm=FALSE,cellnote=dis_distr,
                  notecex=0.1+1/(1.6*log10(nm)),notecol=notecol,na.color="yellow",cexRow = min(1,0.3+1/log10(ns)), 
                  cexCol = min(1.2,0.3+1/log10(nm)),density.info="none",col=color,
                  lhei=c(1,6),lwid=c(1,4),sepwidth=c(0.005,0.005),sepcolor="white",
                  colsep=1:nm,rowsep=1:ns,margins = c(5,6.5))

        dev.off()
    }
}

 
##' generate_angular_distance_trees
##' 
##' build a tree with the angular distance 
##'
##' @export
generate_angular_distance_trees <- function(angular_dist_w_root, sel_normal_sample, colorfile, sample_name_change, ad_tree_w_root_dir, markerset, usedsamps, ns, power, new_sample_names, meanlen_diff_normal, bootstrap, bscut, bsreps, subject_name, markergroup, replicate, ad_matrix_w_root_dir, tree_type='NJ') {

    phylotree <- nj(as.dist(angular_dist_w_root))
    refind <- which(rownames(angular_dist_w_root)==sel_normal_sample)
    phylotree <- root(phylotree,outgroup=refind,resolve.root=TRUE)

    if (bootstrap==FALSE) {
        colors <- colorfile[pmatch(phylotree$tip.label,colorfile[,1]),2]
        ncolor <- paste0("-",colnames(colorfile)[2])
        if (sample_name_change) {
            phylotree$tip.label <- new_sample_names[pmatch(phylotree$tip.label,new_sample_names[,1]),2]

            pdf(file.path(ad_tree_w_root_dir,paste0(subject_name,"_unrooted_angular_distance_tree_newnames_",markergroup,"_",replicate,"_",tree_type,ncolor,".pdf")),width=10,height=10)
            plot.phylo(phylotree,tip.color = colors,font=1,type="unrooted",cex=0.75)
            dev.off() 
        } else {
            pdf(file.path(ad_tree_w_root_dir,paste0(subject_name,"_unrooted_angular_distance_tree_oldnames_",markergroup,"_",replicate,"_",tree_type,ncolor,".pdf")),width=10,height=10)
            plot.phylo(phylotree,tip.color = colors,font=1,type="unrooted",cex=0.75)
            dev.off()
        }

    } else {

        ############
        # bootstrap
        ############

        message('Generating ',bsreps,' bootstrapped angular distance trees/matrices ...')
        bstrees=list()
        bstrees_newnames=list()
        bsmatrices=list()

        for (i in 1:bsreps){
            if(i %% 100 == 0) message('Completed ',i,'/',bsreps,' bootstrap replicates.')
            sample_markers <- sample(markerset,length(markerset),replace=TRUE)
            bs_Z <- Zij(usedsamps,sample_markers,meanlen_diff_normal)
            bs_angular_dist_w_root <- ang_dist_matrix(bs_Z,usedsamps,ns,sample_markers,sel_normal_sample,power)

            ## save the bootstrap angular distance matrix with newnames
            bs_angular_dist_w_root_newnames <- bs_angular_dist_w_root
            newnames <- new_sample_names[pmatch(rownames(bs_angular_dist_w_root_newnames),new_sample_names[,1]),2]
            rownames(bs_angular_dist_w_root_newnames) <- newnames; 
            colnames(bs_angular_dist_w_root_newnames) <- newnames;
            bsmatrices[[i]] <- bs_angular_dist_w_root_newnames

            ## save the bootstrapped NJ tree
            bst <- nj(as.dist(bs_angular_dist_w_root))
            bstrees[[i]] <- bst

            ## make another version with newnames
            bst_newnames <- nj(as.dist(bs_angular_dist_w_root)) 
            bst_newnames$tip.label <- new_sample_names[pmatch(bst_newnames$tip.label,new_sample_names[,1]),2]
            bstrees_newnames[[i]] <- bst_newnames
        }

        ## save .rds object with the bootstrapped trees with newnames
        class(bstrees_newnames) <- "multiPhylo"
        bstrees_file <- file.path(ad_tree_w_root_dir,paste0(subject_name,"_Bootstrapped_unrooted_angular_distance_tree_newnames_",markergroup,"_",replicate,"_",tree_type,"_bstrees.rds"))
        saveRDS(bstrees_newnames,file=bstrees_file)
        bsmatrices_file <- file.path(ad_matrix_w_root_dir,paste0(subject_name,"_angular_dist_matrix_w_root_",markergroup,"_",replicate,"_newnames_bsmatrices.rds"))
        saveRDS(bsmatrices,file=bsmatrices_file)

        ## plot trees with bootstrap values and old sample names
        class(bstrees) <- "multiPhylo"
        par(mar=c(1.1,1.1,1.1,1.1))

        # add some nice colors to the tip labels
        ncolor <- ""
        for (k in 2:ncol(colorfile)) {
            colors <- colorfile[pmatch(phylotree$tip.label,colorfile[,1]),k]
            ncolor <- paste0("-",colnames(colorfile)[k])
            tree <- plotBS(phylotree,bstrees,type="phylogram",p=1)
            bsvalues <- round(tree$node.label,digits=0)
            tree$node.label <- bsvalues
            plot.phylo(tree,tip.color = colors,font=1,cex=0.75)
            nodelabels(tree$node.label,frame="none",adj = c(1, 1.1),cex=0.75)

            pdf(file.path(ad_tree_w_root_dir,paste0(subject_name,"_Bootstrapped_rooted_angular_distance_tree_oldnames_",markergroup,"_",replicate,"_",tree_type,ncolor,".pdf")),width=10,height=10)
            plot.phylo(tree,tip.color = colors,font=1,cex=0.75)
            nodelabels(tree$node.label,frame="none",adj = c(1, 1.1),cex=.75)
            dev.off()

            pdf(file.path(ad_tree_w_root_dir,paste0(subject_name,"_Bootstrapped_unrooted_angular_distance_tree_oldnames_",markergroup,"_",replicate,"_",tree_type,ncolor,".pdf")),width=10,height=10)
            plot.phylo(tree,tip.color = colors,font=1,type="unrooted",cex=0.75)
            nodelabels(tree$node.label,frame="none",adj = c(1, 1.1),cex=.6)
            dev.off()

            pdf(file.path(ad_tree_w_root_dir,paste0(subject_name,"_unrooted_angular_distance_tree_oldnames_",markergroup,"_",replicate,"_",tree_type,ncolor,".pdf")),width=10,height=10)
            plot.phylo(tree,tip.color = colors,font=1,type="unrooted",cex=0.75)
            dev.off()      

            # new sample names
            if (sample_name_change) {
                tree$tip.label <- new_sample_names[pmatch(tree$tip.label,new_sample_names[,1]),2]

                pdf(file.path(ad_tree_w_root_dir,paste0(subject_name,"_Bootstrapped_rooted_angular_distance_tree_newnames_",markergroup,"_",replicate,"_",tree_type,ncolor,".pdf")),width=10,height=10)
                plot.phylo(tree,tip.color = colors,font=1,cex=0.75)
                nodelabels(tree$node.label,frame="none",adj = c(1, 1.1),cex=.75)
                dev.off()

                pdf(file.path(ad_tree_w_root_dir,paste0(subject_name,"_Bootstrapped_unrooted_angular_distance_tree_newnames_",markergroup,"_",replicate,"_",tree_type,ncolor,".pdf")),width=10,height=10)
                plot.phylo(tree,tip.color = colors,font=1,type="unrooted",cex=0.75)
                nodelabels(tree$node.label,frame="none",adj = c(1, 1.1),cex=.6)
                dev.off()

                pdf(file.path(ad_tree_w_root_dir,paste0(subject_name,"_unrooted_angular_distance_tree_newnames_",markergroup,"_",replicate,"_",tree_type,ncolor,".pdf")),width=10,height=10)
                plot.phylo(tree,tip.color = colors,font=1,type="unrooted",cex=0.75)
                dev.off() 
            }
        }
    }
}



