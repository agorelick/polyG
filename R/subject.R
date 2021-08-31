##' subject
##' @export
subject <- function(subject_name,pdir,output_dir,sample_exclusion_th,max_replicate_th,sel_normal_sample,all_normal_samples) { 
    ## - For a given subject, we load the sample-names from the sampleName_file, and read the marker files for the subject 
    ## - then for each marker, we choose a representative replicate.
    ## - Using only one representative replicates for each marker, save the resulting (normalized) data.
    ## - Generate summary stats and plots about the data using the representative replicates for this sample.

    ## define the path to the various poly-G loci data files for this given subject
    data_path <- file.path(pdir,'data',paste0(subject_name,'-Data'))
    data_files <- list.files(data_path)

    ## directory for additional information
    add_info_dir <- file.path(output_dir,"additional_info")
    dir.create(add_info_dir,showWarnings=FALSE )

    ## samples for a subject. Load the sampleName_file if it was specified. Otherwise check for the default file name
    sampleName_file <- file.path("sample_names",paste0(subject_name,"_SampleNames.txt"))
    error_if_file_not_exist(sampleName_file)
    all_sample_names <- sort((read.table(sampleName_file,stringsAsFactors = F,header=F))[,1])
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
    overlap_score <- list()
    df_data_repre_repli <- list() # get marker files with only representative replicates

    i <- 1
    for (marker_file in data_files) {
        df_data <- read.table(file.path(data_path,marker_file),sep="\t",header=T,stringsAsFactors=F,check.names=FALSE)
        #  colnames(df_data) <- gsub(x = colnames(df_data), pattern = "\\.", replacement = "-") 

        #** exclude markers that have only 1 peak
        if (nrow(df_data) < 2) {next}

        marker <- unlist(strsplit(marker_file,".txt"))[1]
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

            ## get the sample name (dropping the '_N' at the end for replicates)
            sample <- paste(head(strsplit(column_name,'_')[[1]],-1),collapse='_')

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

        # create J-S distance matrix
        marker_dm_filename <- paste0(marker,"_distancematrix_df_",max_replicate_th,".tsv")
        raw_dist[[marker]] <- create_JSM_matrix(norm_df,file.path(output_dir,marker_dm_filename))

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

        # get sum of positive components of probability distribution difference between a given sample and normal (reference)
        prob_distr <- norm_df_repre_repli[[marker]]
        overlap_score[[marker]] <- get_overlap_scores(prob_distr,all_sample_names,sel_normal_sample)

        # write all_missing_samples for the marker
        all_missing_samples[[marker]] <- outlist[[marker]]$all_missing_samples
        if (length(all_missing_samples[[marker]]) > 0) {
            write.table(all_missing_samples[[marker]],file=file.path(add_info_dir,paste0(marker,"_all_missing_samples.txt")),sep="\t",quote=F,col.names=F)
        }
        # write samples not in a given marker file
        mis_samples <- outlist[[marker]]$samples_not_in_marker_files
        if (length(mis_samples) > 0) {
            write.table(mis_samples,file=file.path(add_info_dir,paste0(marker,"_samples_not_in_marker_files.txt")),sep="\t",quote=F,col.names=F)
        }
        # write samples without a representative replicate for a given marker
        no_repre_repli_samples <- outlist[[marker]]$samples_wo_repre_repli
        if (length(no_repre_repli_samples) > 0) {
            write.table(no_repre_repli_samples,file=file.path(add_info_dir,paste0(marker,"_samples_wo_repre_repli.txt")),sep="\t",quote=F,col.names=F)
        }

    }
    marker_names <- sort(marker_names)

    # create file with overview about selected replicates and their distance to the normal
    file_name <- paste0(subject_name,"_markers_df_dist_to_normal_",max_replicate_th,".tsv")
    create_overview_table_dist_to_normal(file.path(output_dir,file_name),repre_replicates,raw_dist,all_sample_names,marker_names,sel_normal_sample)

    # create file with JS distance distributions to normal sample
    dis_distr_fp <- file.path(output_dir,paste0(subject_name,"_values_df_",max_replicate_th,".pdf"))
    distance_distribution_plot(dis_distr_fp,repre_replicates,track_lengths,raw_dist,all_sample_names,marker_names,sel_normal_sample,max_replicate_th)

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
    write.table(excluded_samples,file=file.path(add_info_dir,paste0(subject_name,"_excluded_samples_",sample_exclusion_th,".txt")),sep="\t",quote=F,col.names=F)

    # generate distance matrix per marker for the selected replicates (excluded samples not included)
    ma_repdm_filesuff <- paste0("representativedistancematrix_df_",max_replicate_th,".tsv")
    write_distance_matrices(output_dir,ma_repdm_filesuff,sample_names,marker_names,raw_dist,repre_replicates)

    # remove excluded samples from all_missing_samples[[marker]] => missing_sample_ex[[marker]]
    missing_samples_ex <- list()
    for (marker in marker_names) {
        miss <- all_missing_samples[[marker]] 
        miss <- miss[!(miss %in% excluded_samples)]
        missing_samples_ex[[marker]] <- miss  
        if (length(miss) > 0) {
            write.table(missing_samples_ex[[marker]],file=file.path(add_info_dir,paste0(marker,"_missing_samples_exclusion_th_",sample_exclusion_th,".txt")),sep="\t",quote=F,col.names=F)
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
        write(paste0(markers,": ",final_marker_names[[markers]]),file=marker_name_file,append=T)
    }

    # plot overlap scores for all markers per sample
    os_filepath <- file.path(output_dir,paste0(subject_name,"_overlap_scores.txt"))
    osplot_filepath <- file.path(output_dir,paste0(subject_name,"_overlap_scores.pdf"))
    plot_overlap_scores(overlap_score,all_sample_names,marker_names,os_filepath,osplot_filepath)


    ## collect output data to use in sel_rep_replicates_create_dist_matrix.R
    list(nmarkers=nmarkers, final_marker_names=final_marker_names, marker_list=marker_list, marker_names=marker_names, sample_replicates=sample_replicates, repre_replicates=repre_replicates, track_lengths=track_lengths, raw_dist=raw_dist, all_missing_samples=all_missing_samples, outlist=outlist, norm_df_repre_repli=norm_df_repre_repli, overlap_score=overlap_score, df_data_repre_repli=df_data_repre_repli, sample_names=sample_names,add_info_dir=add_info_dir)


}
