ViewAndScoreFull <-
    function(filename,initialExpansion=500,refexpansion=400,mmRate,indelRate,
             bamFilePath='',allowedMM=6,gapOpeningArg = -4,readLength,
             filtsings=TRUE,findSplitReads=FALSE,dedup=TRUE,eventsToEval=Inf,MMsplits=15,pairlimit=2e3,
             gapExtensionArg = -1,substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],
             rngsAlign=GRanges(),build='hg19',bsbuildprefix="BSgenome.Hsapiens.UCSC.",validChr=c(1:22,'X','Y','M'),verbose=FALSE){


        alignedall <- list()
        ## replaced, make sure 'chr' precedes chr #, make sure name of bam MultipleAlnsort.bam in input
        if(is.null(nrow(filename))){
            events <- as.matrix(read.delim(file=filename,header=TRUE))
        }else{
            events <- filename
        }
        events[,'Sample'] <- file.path(bamFilePath,events[,'Sample'])
        if('SplitsSample' %in% colnames(events)) events[,'SplitsSample'] <- file.path(bamFilePath,events[,'SplitsSample'])
        events[,'Chr1'] <- gsub('([Cc][Hh][Rr]){2}','chr',paste('chr',gsub(' ','',events[,'Chr1']),sep=''))
        events[,'Chr2'] <- gsub('([Cc][Hh][Rr]){2}','chr',paste('chr',gsub(' ','',events[,'Chr2']),sep=''))

        if(eventsToEval[1] != Inf){
            events<-events[eventsToEval,]

            if(is.null(ncol(events))) events<-t(events)
        }


        alignedall <- alignViewFull(events,filtsings=filtsings,findSplitReads=findSplitReads,dedup=dedup,initialExpansion=initialExpansion,refexpansion=refexpansion,indelRate=indelRate,mmRate=mmRate,readLength=readLength,allowedMM=allowedMM,gapOpeningArg = gapOpeningArg, gapExtensionArg = gapExtensionArg,substitutionMat=substitutionMat,build=build,rngsAlign=rngsAlign,MMsplits=MMsplits,bsbuildprefix=bsbuildprefix,verbose=verbose)


        if(nrow(events) > 1){
        while(is.list(alignedall) & length(alignedall) != nrow(events)){
            formatted<-unlist(lapply((lapply(alignedall,function(x){c(length(x),length(x[[1]]))})),function(x){all(x==c(3,4))}))
         alignedall <-   c(do.call(c,alignedall[!formatted]),alignedall[formatted])
         if(verbose) print('simplified list object from foreach loop')
        }
    }

        likelihoodScores <- vector()
        forplot <- list()
        if(nrow(events)==1){
            newscore <- Inf; jj <- 1
            cc <- cbind(alignedall[[1]][[3]],alignedall[[2]][[3]],alignedall[[3]][[3]])
            likelis <- alignedall[[1]][[3]]/
                apply(cbind(alignedall[[2]][[3]],alignedall[[3]][[3]]),1,max)
            while(newscore %in% c(Inf,NaN)){
                newscore <- sum(log(likelis,10))
                likelis <- likelis[-which.max(likelis)]
            }
            likelihoodScores[jj] <- newscore
            forplot[[1]] <-list()
            forplot[[1]] <-list(alignedall[[1]][[2]],alignedall[[2]][[2]],alignedall[[3]][[2]])
        }else{
    likelihoodScores=foreach(jj=1:length(alignedall),.combine='c',.multicombine=TRUE) %dopar% {

##            for(jj in 1:length(alignedall)){
                newscore <- Inf
                ## if reads <10^-10 max across all loci are removed, roc is a bit better
                cc <- cbind(alignedall[[jj]][[1]][[3]],alignedall[[jj]][[2]][[3]],alignedall[[jj]][[3]][[3]])
                ##                rem=which(apply(cc,1,max)<10^-6)
                likelis <- alignedall[[jj]][[1]][[3]]/
                    apply(cbind(alignedall[[jj]][[2]][[3]],alignedall[[jj]][[3]][[3]]),1,max)
                ##                if(length(rem)>0) likelis=likelis[-rem]
                while(newscore %in% c(Inf,NaN)){
                    newscore <- sum(log(likelis,10))
                    likelis <- likelis[-which.max(likelis)]
                }
##                likelihoodScores[jj] <- newscore

                newscore
            }
    forplot=lapply(alignedall,function(x){list(x[[1]][[2]],x[[2]][[2]],x[[3]][[2]])})
}

        return(list(forplot=forplot,score=likelihoodScores))
    }
