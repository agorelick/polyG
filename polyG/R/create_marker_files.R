##' create_marker_files
##' @export
create_marker_files <- function(sample_table_file, peak_table_file, markers_file, mixes_file, sizing_artifacts_file=NA, pdir='.', thres=20, intensityfilter=0.1, genemapper=F, overwrite=F) {
   
    ## define the data/ directory expected to be generated. 
    data_dir <- file.path(pdir,"data")
    create_dir_check_if_already_exists(data_dir,overwrite)

    ## define create_marker_file_info directory to store information for troubleshooting
    marker_file_info_dir <- file.path(pdir,"create_marker_file_info")
    create_dir_check_if_already_exists(marker_file_info_dir,overwrite)

    ## Check if the expected input files exist. If not, prompt user to fix and throw error.
    error_if_file_not_exist(sample_table_file)
    error_if_file_not_exist(peak_table_file)
    error_if_file_not_exist(markers_file)
    error_if_file_not_exist(mixes_file)
    if(!is.na(sizing_artifacts_file)) error_if_file_not_exist(sizing_artifacts_file)

    # read all required files 
	samp <- read.table(sample_table_file,sep=",",header=T,stringsAsFactors=F)
	peak <- read.table(peak_table_file,sep=",",header=T,stringsAsFactors=F)
	info <- suppressWarnings(read.table(markers_file,sep="\t",header=T,stringsAsFactors=F))
	mixes <- suppressWarnings(read.table(mixes_file,sep="\t",header=F,stringsAsFactors=F))

    # remove samples with failed sizing
    failedsize <- samp$Sample.File.Name[which(is.na(samp$SQ))]
    if (length(which(is.na(samp$SQ)))!=0) samp <- samp[-which(is.na(samp$SQ)),]

    # sizing artifacts file
    if (!is.na(sizing_artifacts_file)) {
    	sizingartifacts <- T
    	sizearti <- suppressWarnings(read.table(sizing_artifacts_file,sep="\t",header=F,stringsAsFactors=F))
	} else {
		sizingartifacts <- F
	}
	
    # setup list of failed PCR samples
    failedpcr <- list()

    # DEFINE SAMPLES IN PEAK FILE, i.e. for each sample, determine which rows in peak file correspond to its beginning and end
    samples <- peak$Sample.File.Name[1]
    begin <- 1
    end <- 0
    t <- peak$Sample.File.Name[1]

    rec.samples <- for (i in 2:dim(peak)[1]){
        if (peak$Sample.File.Name[i]!=t) {
            end <- c(end,i-1)
            begin <- c(begin,i)
            t <- peak$Sample.File.Name[i]
            samples <- c(samples,t)
    }}

    end <- c(end[-1],dim(peak)[1])

    if (length(which(samples==samp$Sample.File.Name)) != dim(samp)[1]){
        stop("SAMPLES DON'T MATCH UP")
    }

    # extract sample and marker information
    nmarkers <- (dim(info)[1])/2
    ind <- seq(from=1,to=dim(info)[1],by=2)
    markerenu <- list()
    reslist <- list()

    ### k loops through markers
    ### j loops through patients

    # for each marker
    for (k in 1:nmarkers){

        marker <- strsplit(info$X[ind[k]],split="\\.")[[1]][1]
        markerenu[[k]] <- marker
        mix <- mixes[pmatch(marker,mixes[,1]),2]
        dye <- mixes[pmatch(marker,mixes[,1]),3]

        # select names and peak file positions for samples containing that marker
        if(genemapper==F) {
            select <- grep(mix,samp$Sample.Name) 
        } else {
            select <- grep(paste(mix,"\\.",sep=""),samp$Sample.Name)
        }
        ssamples <- samp$Sample.Name[select]
        sbegin <- begin[select]
        send <- end[select]

        patientres <- list()

        #for each patient
        for (j in 2:dim(info)[2]){

            patient <- colnames(info)[j]
            region <- info[c(ind[k],ind[k]+1),j]

            # if a particular marker/allele is not present in a patient, skip
            if (is.na(region[[1]])){
                patientres[[j-1]] <- NA
                next
            }

            # select samples belonging to one patient
            select2 <- grep(patient,ssamples)
            ssamples2 <- ssamples[select2]
            sbegin2 <- sbegin[select2]
            send2 <- send[select2]

            fsizes <- list()
            fheights <- list()
            fareas <- list()

            ############### for each sample
            for (i in 1:length(ssamples2)){

                d <- peak[sbegin2[i]:send2[i],]
                d1 <- d[grep(dye,d$Dye...Sample.Peak),]
                d2 <- d1[which(d1$Size >= region[[1]] & d1$Size <=region[[2]]),]

                # record peaks, make a filler if the marker has not amplified in that sample
                if (length(d2$Size!=0)){

                    fsizes[[i]] <- d2$Size
                    fheights[[i]] <- d2$Height
                    if(genemapper==F) {
                        fareas[[i]] <- d2$Area.Data.Point.
                    } else {
                        fareas[[i]] <- d2$Area
                    }
                } else{
                    fsizes[i] <- NA
                    fheights[i] <- NA
                    fareas[i] <- NA
                }
            }

            ################ end for each sample

            # now get rid of sizing artifacts

            if (sizingartifacts==TRUE){
                samatch <- grep(paste("^",marker,"$",sep=""),sizearti[,1])

                if (length(samatch>0)){
                    if(genemapper==F) {
                        se <- pmatch(paste(mix,"_",sizearti[samatch,2],sep=""),ssamples2)
                    } else {
                        se <- pmatch(paste(mix,sizearti[samatch,2],sep="."),ssamples2)
                    }
                    sampremove <- se[!is.na(se)]

                    if (length(which(!is.na(sampremove)))>0){

                        ssamples2 <- ssamples2[-sampremove]
                        fsizes <- fsizes[-sampremove]
                        fheights <- fheights[-sampremove]
                        fareas <- fareas[-sampremove]
                    }
                }
            }

            defvector <- fsizes[[which.max(lapply(fsizes,length))]]

            # extend defvector as needed to maximum/minimum value in fsizes
            defvector <- c(defvector,seq(from=tail(defvector,n=1),to=(max(unlist(fsizes),na.rm=T)+0.5),by=1)[-1])
            defvector <- c(rev(seq(from=head(defvector,n=1),to=(min(unlist(fsizes),na.rm=T)-0.5),by=-1)[-1]),defvector)

            sizemat <- array(0,dim=c(length(defvector),length(ssamples2)))
            heightmat <- array(0,dim=c(length(defvector),length(ssamples2)))
            areamat <- array(0,dim=c(length(defvector),length(ssamples2)))

            colnames(sizemat) <-ssamples2
            colnames(heightmat) <-ssamples2
            colnames(areamat) <-ssamples2


            #standardize sizes
            nmatch <- function(y,def) which.min(abs(def-y))

            for (y in 1:length(ssamples2)){

                if (length(fsizes[[y]])==1){
                    if (is.na(fsizes[[y]])) next
                }

                ind4 <- sapply(fsizes[[y]],nmatch,def=defvector)
                sizemat[ind4,y] <-fsizes[[y]]
                heightmat[ind4,y] <-fheights[[y]]
                areamat[ind4,y] <-fareas[[y]]

            }

            #REMOVE FAILED PCRs
            failed <- which(apply(heightmat,2,mean)<thres)
            if(length(failed)>0){
                sizemat <- sizemat[,-failed]
                heightmat <- heightmat[,-failed]
                areamat <- areamat[,-failed]
            }


            # filter out PCRs with low intensity
            cutoff <- intensityfilter*mean(apply(heightmat,2,mean))
            exclude <- which(apply(heightmat,2,mean)<cutoff)
            if (length(exclude)>0){
                heightmat <- heightmat[,-exclude]
            }


            #*** create data directory to store marker files for downstream analysis
            patient_data_dir <- file.path(data_dir,paste0(patient,"-Data"))
            dir.create(patient_data_dir,recursive=T,showWarnings=F)

            if(genemapper==F) {
                write.table(heightmat,file.path(patient_data_dir,paste(marker,"_",patient,"_",intensityfilter,".txt",sep="")),sep="\t",quote=F,row.names=F,col.names=as.data.frame(strsplit(colnames(heightmat),paste0(mix,"_")),stringsAsFactors=F)[2,])
            } else {
                write.table(heightmat,file.path(patient_data_dir,paste(marker,"_",patient,"_",intensityfilter,".txt",sep="")),sep="\t",quote=F,row.names=F,col.names=as.data.frame(strsplit(colnames(heightmat),"\\."),stringsAsFactors=F)[2,])
            }

            ## save information for troubleshooting
            write.table(data.frame(Sample=ssamples2,Length=unlist(lapply(fsizes,length))),file=file.path(marker_file_info_dir,paste(marker,"_",patient,"_",intensityfilter,sep="","_peak_length.txt")),sep="\t",quote=F,row.names=F)

        }
        reslist[[k]] <- patientres
        set <- 1

    }
}
