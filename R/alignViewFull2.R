alignViewFull <-
function(events,filtsings=TRUE,findSplitReads=FALSE,dedup=TRUE,allowedMM=6,initialExpansion=0,refexpansion=400,
         indelRate,mmRate,readLength,pairlimit=2e3,gapOpeningArg = -4, gapExtensionArg = -1,
         substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],build='hg19',
         bsbuildprefix="BSgenome.Hsapiens.UCSC.",
         MMsplits=15,rngsAlign=GRanges(),validChr=c(1:22,'X','Y','M'),findSplitReadsOnly=FALSE,verbose=FALSE){


    do.call(library, list(paste(bsbuildprefix, build, sep='')))
    ##pkgname <- paste("BSgenome.Hsapiens.UCSC.", build, sep="")
    genomeName <-gsub("\\.UCSC.+$","",gsub("BSgenome\\.","",bsbuildprefix))

    ii <- NULL
    if(length(gapExtensionArg)==1) gapExtensionArg=rep(gapExtensionArg,nrow(events))
    if(length(gapOpeningArg)==1) gapOpeningArg=rep(gapOpeningArg,nrow(events))
    if(length(allowedMM)==1) allowedMM=rep(allowedMM,nrow(events))
    if(length(initialExpansion)==1) initialExpansion=rep(initialExpansion,nrow(events))
    if(length(refexpansion)==1) refexpansion=rep(refexpansion,nrow(events))
    if(length(rngsAlign)==1 & is.list(rngsAlign)) rngsAlign=rep(list(unlist(rngsAlign)),nrow(events))
    if(is.matrix(substitutionMat)){
        ss=list()
        for(ii in 1:nrow(events)) ss[[ii]]=substitutionMat
        substitutionMat=ss
    }



    alignedall=foreach(ii=1:nrow(events),.combine='list',.multicombine=TRUE) %dopar% {


        bamFile=events[ii,'Sample']
        if('SplitsSample' %in% colnames(events)) bamFileSplits=events[ii,'SplitsSample'] else bamFileSplits <- NULL

        header<-unlist(scanBamHeader(bamFile))
        if(length(grep(paste('SN',events[ii,'Chr1'],sep=':'),header))>0) chrsub <- ''
        else chrsub <- 'chr'
            rngs=GRanges(seqnames=gsub(chrsub,'',c(events[ii,'Chr1'],events[ii,'Chr2'])),
            ranges=IRanges(start=as.numeric(c(events[ii,'Start1'],events[ii,'Start2'])),
            end=as.numeric(c(events[ii,'End1'],events[ii,'End2']))))

        rngs=rngs+initialExpansion[ii]
        if(length(rngsAlign)>0){
            rngsAlign<-rngsAlign[[ii]]+initialExpansion[ii]

        }else{
            rngsAlign<-rngs
        }
        secondrngsAlign=rngsAlign[1]+refexpansion[ii]
        thirdrngsAlign=rngsAlign[2]+refexpansion[ii]


        if(verbose) print(paste('Working on event',ii,'of',nrow(events)))
#        print(c(rngs))
#        print(c(rngsAlign,secondrngsAlign,thirdrngsAlign))
        aligned1=mainAlignViewFull(bamFile,rngs,rngsAlign,filtSings=filtsings,findSplitReads=findSplitReads,dedup=dedup,filterbyMM=TRUE,MM=allowedMM[ii],
        indelRate=indelRate,mmRate=mmRate,readLength=readLength,pairlimit=pairlimit,gapOpeningArg = gapOpeningArg[ii], gapExtensionArg = gapExtensionArg[ii],substitutionMat=substitutionMat[[ii]],bamFileSplits=bamFileSplits,MMsplits=MMsplits,didSplits=FALSE,
        genomeName=genomeName,findSplitReadsOnly=findSplitReadsOnly,doingContig=FALSE,verbose=verbose)
        bamnames=aligned1[[1]]

        if(verbose) print(paste('primary alignment for event',ii,'done'))


        ## added for findSplitReadsOnly
        if(findSplitReadsOnly & length(bamnames)==0){
            alignedcontigs<-list()
            alignedcontigs[[1]]<-alignedcontigs[[2]]<-aligned1
        }else{
            contigrngsAlign<-c(secondrngsAlign,thirdrngsAlign)
        alignedcontigs=mainAlignViewFull(bamFile,rngs,contigrngsAlign,filtSings=filtsings,findSplitReads=findSplitReads,dedup=dedup,filterbyMM=FALSE,
            filterbyname=TRUE,filternames=bamnames,indelRate=indelRate,mmRate=mmRate,readLength=readLength,pairlimit=pairlimit,
            gapOpeningArg = gapOpeningArg[ii],gapExtensionArg = gapExtensionArg[ii],substitutionMat=substitutionMat[[ii]],bamFileSplits=bamFileSplits,
        didSplits=aligned1[[4]],MMsplits=MMsplits,genomeName=genomeName,findSplitReadsOnly=findSplitReadsOnly,doingContig=TRUE,verbose=verbose)
        if(verbose) print(paste('secondary alignment (1 of 2) for event',ii,'done'))

        #aligned3=mainAlignViewFull(bamFile,rngs,thirdrngsAlign,filtSings=filtsings,findSplitReads=findSplitReads,dedup=dedup,filterbyMM=FALSE,filterbyname=TRUE,
        #filternames=bamnames,indelRate=indelRate,mmRate=mmRate,readLength=readLength,pairlimit=pairlimit,gapOpeningArg = gapOpeningArg[ii],
        #gapExtensionArg = gapExtensionArg[ii],substitutionMat=substitutionMat[[ii]],bamFileSplits=bamFileSplits,
        #didSplits=aligned1[[4]],MMsplits=MMsplits,genomeName=genomeName,findSplitReadsOnly=findSplitReadsOnly,verbose=verbose)
        if(verbose) print(paste('secondary alignment (2 of 2) for event',ii,'done'))
    }
        aligned123=list(aligned1,alignedcontigs[[1]],alignedcontigs[[2]])

        return(aligned123)

    }

    return(alignedall)
}
