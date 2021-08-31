##' select_representative_replicates
##' @export
select_representative_replicates <- function(marker,raw_dist,intensity,sample_replicates,repre_replicates,max_replicate_th,all_sample_names,outlist) {
    ### Select most representative replicates. First check if all of them agree, then select the pair
    ### which minimizes the provided distance and then select the replicate that maximizes the intensity

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
                    # all 3 replicates belong to at least one of the pairs with the min(distance)
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
                    # all 4 replicates belong to at least one of the pairs with min(distance)
                    repre_idx <- which(intensity[repli]==max(intensity[repli]))[1]
                    repre_replicates[[marker]][[sample]] <- repli[repre_idx]
                    print(paste0("For ",marker,", ",sample," has more than three pairs of replicates with the minimal distance. Check sizing window."))

                } else if (length(min_dist_idx) == 3) {
                    # out of C(6,3)=20 combinations of 3 distances, all but 4 involves all 4 replicates. The 4 exceptions are distance[c(1,2,5)] - involving replicates (1,2,3); distance[c(1,4,6)] - involving replicates (1,2,4); distance[c(2,3,6)] - involving replicates (2,3,4); distance[3,4,5] - involving replicates (3,4,1)

                    # consider the 4 exceptions first
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
                    # There are C(6,2)=15 combinations of 2 distances. They fall into 5 groups: distance[c(1,2),c(1,5),c(2,5)] involves replicates (1,2,3); distance[c(1,4),c(1,6),c(4,6)] involves replicates (1,2,4); distance[c(2,3),c(2,6),c(3,6)] involves replicates (2,3,4); distance[c(3,4),c(3,5),c(4,5)] involves replicates (1,3,4); distance[c(1,3),c(2,4),c(5,6)] involves replicates (1,2,3,4)

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
                        # the remaining 3 cases involves all 4 replicates
                        repre_idx <- which(intensity[repli]==max(intensity[repli]))[1]
                        repre_replicates[[marker]][[sample]] <- repli[repre_idx]
                        print(paste0("For ",marker,", ",sample," has two pairs of replicates with the minimal distance. Check sizing window."))
                    }

                } else {
                    # length(min_dist_idx)==1
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
