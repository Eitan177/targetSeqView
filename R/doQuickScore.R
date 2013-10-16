doQuickScore <-
    function(filename,
             readLength,bamFilePath='',gapOpeningArg = -4,
             dedup=TRUE,eventsToEval=Inf,mmRate,indelRate,pairlimit=2e3,
             gapExtensionArg = -1,substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],
             rngsAlign=GRanges(),build='hg19',bsbuildprefix="BSgenome.Hsapiens.UCSC.",conservativeContigAlign=FALSE,validChr=c(1:22,'X','Y','M'),verbose=FALSE){

        alignedall <- list()

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

        alignedall <- alignQuick(events,dedup=dedup,indelRate=indelRate,mmRate=mmRate,readLength=readLength,gapOpeningArg = gapOpeningArg,pairlimit=pairlimit, gapExtensionArg = gapExtensionArg,substitutionMat=substitutionMat,build=build,bsbuildprefix=bsbuildprefix,rngsAlign=rngsAlign,conservativeContigAlign=conservativeContigAlign,validChr=validChr,verbose=verbose)

        while(is.list(alignedall) & length(alignedall) != nrow(events)){
         alignedall <-   c(do.call(c,alignedall[unlist(lapply(alignedall,is.list))]),alignedall[!unlist(lapply(alignedall,is.list))])
         if(verbose) print('simplified list object from foreach loop')
        }

        likelihoodScores <- vector()
        if(nrow(events)==1){
            newscore <- Inf; jj <- 1

            likelis <- alignedall[,'pvalSV']/alignedall[,'pvalContig']
            while(newscore %in% c(Inf,NaN)){
                newscore <- sum(log(likelis,10))
                likelis <- likelis[-which.max(likelis)]
            }
            likelihoodScores[jj] <- newscore
        }else{
            likelihoodScores=foreach(jj=1:length(alignedall),.combine='c',.multicombine=TRUE) %dopar% {
                newscore <- Inf
                likelis <- alignedall[[jj]][,'pvalSV']/alignedall[[jj]][,'pvalContig']
                while(newscore %in% c(Inf,NaN)){
                    newscore <- sum(log(likelis,10))
                    likelis <- likelis[-which.max(likelis)]
                }
                ##print(newscore)
                newscore
            }
        }
        if(eventsToEval[1] != Inf){
            names(likelihoodScores) <- eventsToEval
        }else{
            names(likelihoodScores) <- seq_along(likelihoodScores)
        }
        return(likelihoodScores)
    }
