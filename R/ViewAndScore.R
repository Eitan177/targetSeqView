ViewAndScore <-
    function(filename,initialExpansion=500,refexpansion=400,estimateIndelRate=TRUE, indelRate=0.005,estimateMmRate=TRUE,mmRate=0.01,
             getReadLength=TRUE,readLength=100,scoreOnly=FALSE,normalBam='',bamFilePath='',allowedMM=6,gapOpeningArg = -4, gapExtensionArg = -1,
             substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],build='hg19',verbose=FALSE){

        alignedall <- list()
        ## replaced, make sure 'chr' precedes chr #, make sure name of bam MultipleAlnsort.bam in input
        events <- as.matrix(read.delim(file=filename,header=TRUE))
        events[,'Sample'] <- file.path(bamFilePath,events[,'Sample'])
        events[,'Chr1'] <- gsub('chrchr','chr',tolower(paste('chr',gsub(' ','',events[,'Chr1']),sep='')))
        events[,'Chr2'] <- gsub('chrchr','chr',tolower(paste('chr',gsub(' ','',events[,'Chr2']),sep='')))
        if(estimateIndelRate || estimateMmRate || getReadLength){
            norms <- normalAlign(bamFile=normalBam,build=build)
            rL <- max(width((DNAStringSet(gsub('-','',as.character(pattern(norms)))))))
            readLength <- ifelse(getReadLength,rL,readLength)
            indls <- nindel(norms)
            indelRateEstimate <- mean(indls@insertion[,2]+indls@deletion[,1])/readLength
            indelRate <- ifelse(estimateIndelRate,indelRateEstimate,indelRate)
            print(indelRate)
            mms <- rep(0,length(norms))
            mms[as.numeric(names(table(mismatchTable(norms)[,1])))] <- table(mismatchTable(norms)[,1])
            mmRateEstimate <- mean(mms)/readLength
            mmRate <- ifelse(estimateMmRate,mmRateEstimate,mmRate)
            print(mmRate)
        }
        alignedall <- alignView(events,scoreOnly=scoreOnly,initialExpansion=initialExpansion,refexpansion=refexpansion,indelRate=indelRate,mmRate=mmRate,readLength=readLength,allowedMM=allowedMM,gapOpeningArg = gapOpeningArg, gapExtensionArg = gapExtensionArg,substitutionMat=substitutionMat,build=build,verbose=verbose)

        likelihoodScores <- vector()
        for(jj in 1:length(alignedall)){
            newscore <- Inf
            ## if reads <10^-10 max across all loci are removed, roc is a bit better
            cc <- cbind(alignedall[[jj]][[1]][[3]],alignedall[[jj]][[2]][[3]],alignedall[[jj]][[3]][[3]])
            ##                rem=which(apply(cc,1,max)<10^-6)
            likelis <- alignedall[[jj]][[1]][[3]]/
                apply(cbind(alignedall[[jj]][[2]][[3]],alignedall[[jj]][[3]][[3]]),1,max)
            ##                if(length(rem)>0) likelis=likelis[-rem]
            while(newscore %in% c(Inf,NaN)){
                newscore <- sum(log(likelis))
                likelis <- likelis[-which.max(likelis)]
            }
            likelihoodScores[jj] <- newscore
        }

        return(list(alignedall,likelihoodScores))
    }
