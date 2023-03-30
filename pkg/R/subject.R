##' subject
##' 
##' Process the raw poly-G marker data for this subject. Results will be used to in angular_distance() function.
##'
##' @export
subject <- function(input_dir, allout_dir, subject_name, data_path, sample_exclusion_th, sel_normal_sample, all_normal_samples, new_sample_names) { 

    ## raw input data file for each marker
    data_files <- dir(file.path(input_dir,"data",data_path), full.names=T)

    ## for output
    output_dir <- file.path(allout_dir,paste0(subject_name,"_R"))
    if (!dir.exists(output_dir)) dir.create(output_dir,recursive=TRUE,showWarnings=FALSE )

    ## directory for additional information
    add_info_dir <- file.path(output_dir,"additional_info")
    dir.create(add_info_dir,showWarnings=FALSE )

    ## samples for a subject
    sampleName_file <- file.path(input_dir, "sample_names",paste0(subject_name,"_SampleNames.txt"))
    all_sample_names <- sort((read.table(sampleName_file,stringsAsFactors = F))[,1])
    nsamples <- length(all_sample_names)

    ## read marker files
    nmarkers <- length(data_files)
    marker_names <- character()
    sample_replicates <- list()
    repre_replicates <- list()
    track_lengths <- list()
    raw_dist <- list()
    all_missing_samples <- list()
    outlist <- list()
    norm_df_repre_repli <- list()
    df_data_repre_repli <- list() # get marker files with only representative replicates

    i <- 1
    for (marker_file in data_files) {
        df_data <- read.table(marker_file,sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)

        #** exclude markers that have only 1 peak
        if (nrow(df_data) < 2) {next}

        marker <- tail(strsplit(marker_file,'[/]')[[1]], 1)
        marker <- strsplit(marker,'[_]')[[1]][1] #gsub('[.]txt','',marker)
        marker_names[i] <- marker
        i <- i + 1

        # normalize data => each column (sample/replicate) sums to 1
        intensity <- apply(df_data,2,function(x){sum(x)})
        #** remove markers which has one or more columns of 0s. They should've been filtered out.
        if (length(which(intensity==0)) > 0) {
            print(paste0("Some columns in ",marker_file," are all 0s! This marker is not used."))
            next
        }
        norm_df <- apply(df_data,2,function(x){x/sum(x)})

        nrows <- dim(norm_df)[1]
        ncols <- dim(norm_df)[2]

        # samples and their replicates in the marker file 
        sample_replicates[[marker]] <- list()
        repre_replicates[[marker]] <- list()

        for (sample in all_sample_names) {
            sample_replicates[[marker]][[sample]] <- character()
            repre_replicates[[marker]][[sample]] <- character()
        }
        for (col_idx in 1:ncols) {
            column_name <- colnames(norm_df)[col_idx]
            #sample <- unlist(strsplit(column_name,"_"))[1]
            # assuming <= 10 replicates which are denoted by a separator and a single alphanumeric character after sample name
            sample <- substr(column_name,1,nchar(column_name)-2)
            # check column names against all_sample_names
            if (!(sample %in% all_sample_names)) {
                cat(paste0("Error: incorrect sample name, ",sample," in ",subject_name," marker ",marker,"\n"))
            }
            sample_replicates[[marker]][[sample]] <- append(sample_replicates[[marker]][[sample]],column_name)
        }

        # calculate polyG length statistics - track_lengths for each marker only defined for sample replicates in that marker file
        track_lengths[[marker]] <- list()
        for (column_name in colnames(norm_df)) {
            lengths <- c()
            for (row_idx in 1:nrows) {
                lengths <- c(lengths,rep(row_idx,df_data[row_idx,column_name]))
            }
            len_mean <- mean(lengths)
            len_median <- median(lengths)

            df_data_col <- df_data[,column_name]
            m <- max(df_data_col)   
            lpeak <- which(df_data_col==m)[1]  
            m2 <- sort(df_data_col,partial=nrows-1)[nrows-1]
            lpeak2 <- which(df_data_col==m2)[1]
            mean3_lowIdx <- 1
            if (lpeak > 2) {
                for (idx in 1:(lpeak-2)) {
                    mean3_lowIdx <- mean3_lowIdx + df_data_col[idx]
                }
            }
            len <- length(lengths)
            mean3_hiIdx <- len
            if ((lpeak+1) < nrows) {
                mean3_hiIdx <- 0
                for (idx in 1:(lpeak+1)) {
                    mean3_hiIdx <- mean3_hiIdx + df_data_col[idx]
                }
            }
            len_mean3 <- mean(lengths[mean3_lowIdx:mean3_hiIdx])
            track_lengths[[marker]][[column_name]] <- list(lmean=len_mean,lmedian=len_median,lpeak=lpeak,lpeak2=lpeak2,lmean3=len_mean3)
        }
        raw_dist[[marker]] <- create_JSM_matrix(norm_df)

        # select representative replicate of each sample
        outlist[[marker]] <- list()
        outlist[[marker]] <- select_representative_replicates(marker,raw_dist,intensity,sample_replicates,repre_replicates,max_replicate_th,all_sample_names,outlist)
        repre_replicates[[marker]] <- outlist[[marker]]$repre_replicates

        # get norm_df_repre_repli from the columns of norm_df that correpsond to the representative replicates for each sample in the marker file
        norm_df_repre_repli[[marker]] <- data.frame()
        all_repre_replicates <- unlist(repre_replicates[[marker]])
        sel <- which(colnames(norm_df) %in% all_repre_replicates)
        norm_df_repre_repli[[marker]] <- norm_df[,sel]
        colnames(norm_df_repre_repli[[marker]]) <- sapply(colnames(norm_df_repre_repli[[marker]]),function(x) substr(x,1,nchar(x)-2))

        # write normalized_data (representative replicates only) and normalized_marker files (all replicates) for future reference
        normalized_data_dir <- file.path(output_dir,"normalized_data")
        dir.create(normalized_data_dir,recursive=TRUE,showWarnings=FALSE)
        write.table(norm_df_repre_repli[[marker]],file=file.path(normalized_data_dir,paste0(marker,"_normalized.tsv")),sep="\t",quote=FALSE,row.names=FALSE)

        normalized_marker_files_dir <- file.path(output_dir,"normalized_marker_files")
        dir.create(normalized_marker_files_dir,recursive=TRUE,showWarnings=FALSE)
        write.table(norm_df,file=file.path(normalized_marker_files_dir,paste0(marker,"_normalized_wReplicates.tsv")),sep="\t",quote=FALSE,row.names=FALSE)

        # get df_data_repre_repli from the columns of df_data that correpsond to the representative replicates for each sample in the marker file
        df_data_repre_repli[[marker]] <- data.frame()
        df_data_repre_repli[[marker]] <- df_data[,sel]

        repre_repli_data_dir <- file.path(output_dir,"repre_repli_data")
        dir.create(repre_repli_data_dir,recursive=TRUE,showWarnings=FALSE)
        write.table(df_data_repre_repli[[marker]],file=file.path(repre_repli_data_dir,paste0(marker,"_repre_repli.txt")),sep="\t",quote=FALSE,row.names=FALSE)

        # write all_missing_samples for the marker
        all_missing_samples[[marker]] <- outlist[[marker]]$all_missing_samples
        if (length(all_missing_samples[[marker]]) > 0) {
            write.table(all_missing_samples[[marker]],file=file.path(add_info_dir,paste0(marker,"_all_missing_samples.txt")),sep="\t",quote=FALSE,col.names=F)
        }
        # write samples not in a given marker file
        mis_samples <- outlist[[marker]]$samples_not_in_marker_files
        if (length(mis_samples) > 0) {
            write.table(mis_samples,file=file.path(add_info_dir,paste0(marker,"_samples_not_in_marker_files.txt")),sep="\t",quote=FALSE,col.names=F)
        }
        # write samples without a representative replicate for a given marker
        no_repre_repli_samples <- outlist[[marker]]$samples_wo_repre_repli
        if (length(no_repre_repli_samples) > 0) {
            write.table(no_repre_repli_samples,file=file.path(add_info_dir,paste0(marker,"_samples_wo_repre_repli.txt")),sep="\t",quote=FALSE,col.names=F)
        }

    }
    marker_names <- sort(marker_names)

    # number of unavailable markers for a given sample
    nmismarkers <- list()
    for (sample in all_sample_names) {
        nmismarkers[[sample]] <- 0
    }

    for (marker in marker_names) {
        if (length(all_missing_samples[[marker]]) > 0) {
            for (sample in all_missing_samples[[marker]]) {
                nmismarkers[[sample]] <- nmismarkers[[sample]] + 1
            }
        }
    }

    # exclude samples missing more than sample_exclusion_th of all markers
    # also exclude additional normal samples that miss 1 or more markers so that they don't cause more markers to be discarded
    sample_names <- all_sample_names
    excluded_samples <- character()
    for (sample in all_sample_names) {
        if ((sample %in% setdiff(all_normal_samples,sel_normal_sample)) & (nmismarkers[[sample]] > 0)) {
            excluded_samples <- append(excluded_samples, sample)
            sample_names <- setdiff(sample_names,sample)
            next
        }
        if (nmismarkers[[sample]] > sample_exclusion_th * nmarkers) {
            excluded_samples <- append(excluded_samples, sample)
            sample_names <- setdiff(sample_names,sample)
        }
    }
    sample_names <- sort(sample_names)
    write.table(excluded_samples,file=file.path(add_info_dir,paste0(subject_name,"_excluded_samples_",sample_exclusion_th,".txt")),sep="\t",quote=FALSE,col.names=F)


    # generate distance matrix per marker for the selected replicates (excluded samples not included)
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

            write.table(ma_dist_df,file=file.path(output_dir,'representative_distance_matrices',paste0(marker,"_",fn_suffix)),sep="\t",quote=FALSE,col.names=NA)
        }
    }
    ma_repdm_filesuff <- paste0("representativedistancematrix_df_",max_replicate_th,".tsv")
    dir.create(file.path(output_dir,'representative_distance_matrices'))
    write_distance_matrices(output_dir,ma_repdm_filesuff,sample_names,marker_names,raw_dist,repre_replicates)


    # remove excluded samples from all_missing_samples[[marker]] => missing_sample_ex[[marker]]
    missing_samples_ex <- list()
    for (marker in marker_names) {
        miss <- all_missing_samples[[marker]] 
        miss <- miss[!(miss %in% excluded_samples)]
        missing_samples_ex[[marker]] <- miss  
        if (length(miss) > 0) {
            write.table(missing_samples_ex[[marker]],file=file.path(add_info_dir,paste0(marker,"_missing_samples_surviving_exclusion_th_",sample_exclusion_th,".txt")),sep="\t",quote=FALSE,col.names=FALSE)
        }
    }

    # record all markes and used markers
    marker_name_file <- file.path(add_info_dir,paste0(subject_name,"_final_marker_names.txt"))
    if (file.exists(marker_name_file)) {
        file.remove(marker_name_file)
    }
    final_marker_names <- list()
    marker_list <- c("all_markers","used_markers")
    for (markers in marker_list) {
        final_marker_names[[markers]] <- marker_names
        if (markers == "used_markers") {
            #  remove markers with missing sample representative replicates
            for (marker in marker_names) {
                if (length(missing_samples_ex[[marker]]) > 0) {
                    final_marker_names[[markers]] <- setdiff(final_marker_names[[markers]],marker)
                }
            }
        }
        write(paste0(markers,": ",final_marker_names[[markers]]),file=marker_name_file,append=TRUE)
    }

}




