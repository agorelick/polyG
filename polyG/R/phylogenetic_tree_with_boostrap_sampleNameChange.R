##' phylogenetic_tree_with_boostrap_sampleNameChange
##' @export
phylogenetic_tree_with_boostrap_sampleNameChange <- function(subject_name,pdir,output_dir,sample_exclusion_th,max_replicate_th,sel_normal_sample,all_normal_samples,colorfile_header,KLM,allout_dir,bscut,method) {

    ### 7/19/2020 Add the calculations of the normalized pairwise internal node distances between all pairs of samples
    ### 7/16/2020 New sample names are changed to without space. Need to rewrite the functions.
    ### 7/14/2020 Move rds_klm_onemet function to this script to calculate KLM for inter-lesion mets in all-sample trees
    ###     (rds_klm_mettype function stays in KLM script to calculate KLM for mtes in collapsed trees)
    ###
    ### This script plots both trees with old sample names and those with new sample names
    ### plot old ones first, then plot new ones with tip labels changed
    ###
    ### Required: the final sample names (2nd column of the annotation files) should be
    ###    of the form of "tumor type" followed by number, followed by optional lower case letter
    ###    e.g. Per1, Per10a, Per10b, etc. where Per1 is a single-sample lesion,
    ###         while Per10 is a multi-sample lesion, and Per10a and Per10b denote 
    ###         two samples from the same lesion Per10

    # load required packages
    require(ape)
    require(phangorn)
    require(adephylo)
    
	if(!method %in% c('UPGMA','NJ')) stop('method argument must be either "NJ" or "UPGMA"')

    new_sample_names <- read.table(file.path(pdir,"AnnotationFiles",paste0(subject_name,"_annotationFile.txt")),sep="\t",header=TRUE,stringsAsFactors=FALSE)

    ###############################################
    rds_klm_onemet <- function(phy){
        phy$tip.label <- new_sample_names[pmatch(phy$tip.label,new_sample_names[,1]),2]

        # get sample type of each sample
        st <- sapply(phy$tip.label,function(x) unlist(strsplit(x,"[0-9]+"))[1])
        #st <- sapply(phy$tip.label,function(x) unlist(strsplit(x," "))[1])
        #* change "PerOv" to "Per"
        st[which(st=="PerOv")] <- "Per"
        #* change all primary tumors to paper_primary_tumor_label (e.g., PTa, PTb, PTc to PT)
        st[grep(paper_primary_tumor_label,st)] <- paper_primary_tumor_label

        # select tumor samples
        sel_t <- st!=paper_normal_sample_label
        tt <- st[sel_t]

        # number of tumor samples
        nts <- length(tt)

        # select metastatic samples
        sel_m <- !(st %in% c(paper_primary_tumor_label,paper_normal_sample_label))

        # get the name of metastasis (lesion) for each sample (e.g., Per 1)
        samp_mn <- sapply(phy$tip.label[sel_m],function(x) strsplit(x,"[a-z]*$")[[1]])

        # get the set of names of the metastases
        met_name <- unique(samp_mn)

        # ignore trees in which every metastasis has only one sample
        if (length(met_name) < length(samp_mn)) {
            # get number of samples for each metastasis (met_name): m
            nmn <- integer()
            for (mn in met_name) {
                nmn[mn] <- length(which(samp_mn==mn))
            }

            k <- integer()
            l <- integer()
            m <- integer()

            rootID <- length(phy$tip.label) + 1
            for (mn in met_name) {
                # find klm for each metastasis that has more than one sample
                if (nmn[mn] < 2) next
                m[mn] <- nmn[mn]
                k[mn] <- nts - m[mn]

                l[mn] <- 1
                for (i in rootID:max(phy$edge[,1])) {
                    clade <- extract.clade(phy, i)
                    labels <- clade$tip.label
                    clade_size <- length(labels)

                    # get the size of the largest clade formed by a given met_name
                    labels_m <- sapply(labels,function(x) strsplit(x,"[a-z]*$")[[1]])
                    if (length(which(labels_m==mn))==clade_size) {
                        if (clade_size > l[mn]) {
                            l[mn] = clade_size
                        }
                    }
                }
            }
            df = data.frame(k=k,l=l,m=m)
            return(df)
        } else {
            df = data.frame(c=0)
            return(df)
        }
    }	
    ###############################################

    # define a variable called "folder" which will be the step 2 output folder for that patient
    folder <- file.path(output_dir)
    
    colorfile <- read.table(file.path(pdir,"colorfiles",paste0(subject_name,"_colorfile.txt")),sep="\t",header=colorfile_header,stringsAsFactors=FALSE,fill=TRUE)
    datamat <- grep('representativedistancematrix_df',dir(folder,full.names=T),value=T)
    
    # collect marker information
    matrices <- list()
    matrixcounter <- 0

    for (i in 1:length(datamat)) {

        m <- read.table(datamat[i],sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,na.strings="NA")

        # skip markers that were excluded
        if (sum(is.na(m[,1])) > 0) next

        matrixcounter <- matrixcounter+1
        matrices[[matrixcounter]] <- m

        if (matrixcounter==1){
            allsum <- m
        }else{
            allsum <- allsum+m
        }
    }

    # normalize sum to the number of markers used
    allsum <- allsum/length(matrices)

    #** when sample names have a "-", they are changed to "." in colnames
    colnames(allsum) <- rownames(allsum)
    # where is the normal reference?
    refind <- which(rownames(allsum)==sel_normal_sample)
    #refind <- pmatch(ref,colnames(allsum))

    # the resulting phylogenetic tree    
    phylotree <- nj(as.dist(allsum))
    phylotree <- root(phylotree,outgroup=refind,resolve.root=TRUE)

    if (KLM) {
        klm_onemet <- list()
        #* get klm values for each metastasis with more than one sample
        klm_per_meta <- rds_klm_onemet(phylotree)
        if (klm_per_meta[1,1] > 0) {
            klm_onemet[[subject_name]] <- klm_per_meta
            # the name of metastasis with more than one sample
            metname <- rownames(klm_onemet[[subject_name]])
            KLM_onemet_dir <- file.path(pdir,allout_dir,"KLM_onemet")
            dir.create(KLM_onemet_dir,recursive = TRUE, showWarnings=FALSE )
            write.table(klm_onemet[[subject_name]],file=file.path(KLM_onemet_dir,paste0(subject_name,"_klm_onemet_NJ.txt")),sep="\t",quote=FALSE,col.names=NA)
        }
    }


    #* get normalized pairwise internal node distance 
    nodedist <- as.matrix(distTips(phylotree,tips="all",method="nNodes"))
    normfactor <- phylotree$Nnode
    colnames(nodedist) <- new_sample_names[match(phylotree$tip.label,new_sample_names[,1]),2]
    rownames(nodedist) <- colnames(nodedist)
    norm_nodedist <- nodedist/normfactor
    norm_nodedist_dir <- file.path(pdir,allout_dir,"normalized_internal_node_distance_AllSampleTrees")
    dir.create(norm_nodedist_dir,recursive = TRUE, showWarnings=FALSE )
    write.table(norm_nodedist,file=file.path(norm_nodedist_dir,paste0(subject_name,"_normalized_internal_node_distance_AllSampleTree_NJ.txt")),sep="\t",quote=FALSE,col.names=NA)

    ############
    # bootstrap
    ############

    reps=1000
    bstrees=list()
    bstreesnames=list()

    if (KLM) {
        bsklm_onemet=list()  # data frame with metastasis names as rownames
        bsklm_onemet_mn=list() # data frame to store 1000 bootstrap klm values for each individual metastasis

        if (klm_per_meta[1,1] > 0) {
            for (mn in metname) {
                bsklm_onemet_mn[[mn]] <- c()
            }
        }
    }


    for (i in 1:reps){
        sample <- sample(length(matrices),length(matrices),replace=TRUE)
        bsmat <- matrices[[sample[1]]]
        for (j in 2:length(sample)) bsmat <- bsmat+matrices[[sample[j]]]
        if(method=='UPGMA') {
        	bst <- upgma(as.dist(bsmat))
        } else if(method=='NJ') {
            bst <- nj(as.dist(bsmat))
        }
        bstrees[[i]] <- bst
        bstreesnames[[i]] <-bst

        if (KLM) {
            # do not need this for trees in which every metastasis has only one sample
            if (klm_per_meta[1,1] > 0) {
                #metname <- rownames(klm_per_meta)
                bsklm_onemet[[i]] <- rds_klm_onemet(bstrees[[i]])

                # for each metastasis (lesion) with more than one sample, get a row of klm for each bootstrp repetition
                for (mn in metname) {
                    bsklm_onemet_mn[[mn]] <- rbind(bsklm_onemet_mn[[mn]],as.integer(bsklm_onemet[[i]][mn,]))
                }
            }
        }

    }
    if (KLM) {
        # write bsklm output files  
        if (klm_per_meta[1,1] > 0) {
            for (mn in metname) {
                bsklm_onemet_dir <- file.path(pdir,allout_dir,"bsklm_onemet")
                dir.create(bsklm_onemet_dir,recursive = TRUE, showWarnings=FALSE)
                write.table(bsklm_onemet_mn[[mn]],file=file.path(bsklm_onemet_dir,paste0(subject_name,"_",mn,"_bsklm_onemet_NJ.txt")),sep="\t",quote=FALSE,col.names=c("k","l","m"),row.names=FALSE)
            }
        }
    }

    ## plot trees with bootstrap values and old sample names

    class(bstrees) <- "multiPhylo"
    par(mar=c(1.1,1.1,1.1,1.1))

    # add some nice colors to the tip labels
    ncolor <- ""
    for (k in 2:ncol(colorfile)) {
        colors <- colorfile[pmatch(phylotree$tip.label,colorfile[,1]),k]
        if (colorfile_header) {
            ncolor <- paste0("-",colnames(colorfile)[k])
        }

        tree <- plotBS(phylotree,bstrees,type="phylogram",p=1)
        bsvalues <- round(tree$node.label,digits=0)
        tree$node.label <- bsvalues
        plot.phylo(tree,tip.color = colors,font=1,cex=0.75)
        nodelabels(tree$node.label,frame="none",adj = c(1, 1.1),cex=0.75)

        pdf(paste0("Bootstrapped_Rooted_AllSample_Tree_",tail(strsplit(folder,"\\/")[[1]],n=1),ncolor,"_oldnames_NJ.pdf"),width=9,height=7)
        plot.phylo(tree,tip.color = colors,font=1,cex=0.75)
        nodelabels(tree$node.label,frame="none",adj = c(1, 1.1),cex=0.75)
        dev.off()

        pdf(paste0("Bootstrapped_Unrooted_AllSample_Tree_",tail(strsplit(folder,"\\/")[[1]],n=1),ncolor,"_oldnames_NJ.pdf"),width=9,height=7)
        plot.phylo(tree,tip.color = colors,font=1,type="unrooted",cex=0.75)
        nodelabels(tree$node.label,frame="none",adj = c(1, 1.1),cex=.75)
        dev.off()

        pdf(paste0("Unrooted_AllSample_Tree_with","_",tail(strsplit(folder,"\\/")[[1]],n=1),ncolor,"_oldnames_NJ.pdf"),width=9,height=7)
        plot.phylo(tree,tip.color = colors,font=1,type="unrooted",cex=0.75)
        dev.off()

        options(warn=-1)    
        tree2 <- pruneTree(tree,bscut)
        options(warn=0)

        pdf(paste0("Bootstrapped_Rooted_Pruned_AllSample_Tree_",tail(strsplit(folder,"\\/")[[1]],n=1),"_bscut",bscut,ncolor,"_oldnames_NJ.pdf"),width=9,height=7)
        plot.phylo(tree2,tip.color = colors,font=1,cex=0.75)
        nodelabels(tree2$node.label,frame="none",adj = c(1, 1.1),cex=0.75)
        dev.off()

        ## plot trees with bootstrap values and new sample names

        #** change sample names (colnames and rownames) 
        #new_sample_names <- read.table(file.path(pdir,"AnnotationFiles",paste0(subject_name,"_annotationFile.txt")),sep="\t",header=TRUE,stringsAsFactors=FALSE)
        #ref_new <- new_sample_names[pmatch(ref,new_sample_names[,1]),2]


        #class(bstrees) <- "multiPhylo"
        #par(mar=c(1.1,1.1,1.1,1.1))
        #tree <- plotBS(phylotree,bstrees,type="phylogram",p=1)
        tree$tip.label <- new_sample_names[pmatch(tree$tip.label,new_sample_names[,1]),2]
        #bsvalues <- round(tree$node.label,digits=0)
        #tree$node.label <- bsvalues
        plot.phylo(tree,tip.color = colors,font=1,cex=0.75)
        nodelabels(tree$node.label,frame="none",adj = c(1, 1.1),cex=0.75)

        pdf(paste0("Bootstrapped_Rooted_AllSample_Tree_",tail(strsplit(folder,"\\/")[[1]],n=1),ncolor,"_newnames_NJ.pdf"),width=9,height=7)
        plot.phylo(tree,tip.color = colors,font=1,cex=0.75)
        nodelabels(tree$node.label,frame="none",adj = c(1, 1.1),cex=0.75)
        dev.off()

        pdf(paste0("Bootstrapped_Unrooted_AllSample_Tree_",tail(strsplit(folder,"\\/")[[1]],n=1),ncolor,"_newnames_NJ.pdf"),width=9,height=7)
        plot.phylo(tree,tip.color = colors,font=1,type="unrooted",cex=0.75)
        nodelabels(tree$node.label,frame="none",adj = c(1, 1.1),cex=0.75)
        dev.off()

        pdf(paste0("Unrooted_AllSample_Tree_with","_",tail(strsplit(folder,"\\/")[[1]],n=1),ncolor,"_newnames_NJ.pdf"),width=9,height=7)
        plot.phylo(tree,tip.color = colors,font=1,type="unrooted",cex=0.75)
        dev.off()

        options(warn=-1)    
        tree2 <- pruneTree(tree,bscut)
        options(warn=0)

        pdf(paste0("Bootstrapped_Rooted_Pruned_AllSample_Tree_",tail(strsplit(folder,"\\/")[[1]],n=1),"_bscut",bscut,ncolor,"_newnames_NJ.pdf"),width=9,height=7)
        plot.phylo(tree2,tip.color = colors,font=1,cex=0.75)
        nodelabels(tree2$node.label,frame="none",adj = c(1, 1.1),cex=0.75)
        dev.off()
    }

}

