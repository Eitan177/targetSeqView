alignView <-
<<<<<<< HEAD
function(events,filtsings=TRUE,scoreOnly=FALSE,allowedMM=6,initialExpansion=0,refexpansion=400,
=======
function(events,filtsings=TRUE,findSplitReads=FALSE,dedup=TRUE,scoreOnly=FALSE,allowedMM=6,initialExpansion=0,refexpansion=400,
>>>>>>> fe8fe54aa12fa5f53a9270a3804121d5b15455e9
         indelRate=.005,mmRate=.01,readLength=100,gapOpeningArg = -4, gapExtensionArg = -1,
         substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],build='hg19',verbose=FALSE){

    do.call(library, list(paste("BSgenome.Hsapiens.UCSC.", build, sep='')))
    pkgname <- paste("BSgenome.Hsapiens.UCSC.", build, sep="")
    ii <- NULL
<<<<<<< HEAD
=======
    if(length(gapExtensionArg)==1) gapExtensionArg=rep(gapExtensionArg,nrow(events))
    if(length(gapOpeningArg)==1) gapOpeningArg=rep(gapOpeningArg,nrow(events))
    if(length(allowedMM)==1) allowedMM=rep(allowedMM,nrow(events))
    if(length(initialExpansion)==1) initialExpansion=rep(initialExpansion,nrow(events))
    if(length(refexpansion)==1) refexpansion=rep(refexpansion,nrow(events))
    if(is.matrix(substitutionMat)){
        ss=list()
        for(ii in 1:nrow(events)) ss[[ii]]=substitutionMat
        substitutionMat=ss
    }

>>>>>>> fe8fe54aa12fa5f53a9270a3804121d5b15455e9
    alignedall=foreach(ii=1:nrow(events),.combine='list',.multicombine=TRUE) %dopar% {


        bamFile=events[ii,'Sample']
        rngs=GRanges(seqnames=gsub('chr','',c(events[ii,'Chr1'],events[ii,'Chr2'])),
        ranges=IRanges(start=as.numeric(c(events[ii,'Start1'],events[ii,'Start2'])),
        end=as.numeric(c(events[ii,'End1'],events[ii,'End2']))))
<<<<<<< HEAD
        rngs=rngs+initialExpansion
        secondrngs=rngs[1]+refexpansion
        thirdrngs=rngs[2]+refexpansion

        if(verbose) print(paste('Working on event',ii,'of',nrow(events)))

        aligned1=mainAlignView(bamFile,rngs,rngs,filtSings=filtsings,filterbyMM=TRUE,MM=allowedMM,returnScoreOnly=scoreOnly,
        indelRate=indelRate,mmRate=mmRate,readLength,gapOpeningArg = gapOpeningArg, gapExtensionArg = gapExtensionArg,substitutionMat=substitutionMat)
        bamnames=aligned1[[1]]

        if(verbose) print(paste('primary alignment for event',ii,'done'))
                                        #dev.set(3)
        aligned2=mainAlignView(bamFile,rngs,secondrngs,filtSings=filtsings,filterbyMM=FALSE,filterbyname=TRUE,
        returnScoreOnly=scoreOnly,filternames=bamnames,indelRate=indelRate,mmRate=mmRate,readLength,gapOpeningArg = gapOpeningArg, gapExtensionArg = gapExtensionArg)
        if(verbose) print(paste('secondary alignment (1 of 2) for event',ii,'done'))
                                        #dev.set(4)
        aligned3=mainAlignView(bamFile,rngs,thirdrngs,filtSings=filtsings,filterbyMM=FALSE,filterbyname=TRUE,
        returnScoreOnly=scoreOnly,filternames=bamnames,indelRate=indelRate,mmRate=mmRate,readLength,gapOpeningArg = gapOpeningArg, gapExtensionArg = gapExtensionArg)
        if(verbose) print(paste('secondary alignment (2 of 2) for event',ii,'done'))
                                        #browser()
=======
        rngs=rngs+initialExpansion[ii]
        secondrngs=rngs[1]+refexpansion[ii]
        thirdrngs=rngs[2]+refexpansion[ii]

        if(verbose) print(paste('Working on event',ii,'of',nrow(events)))

        aligned1=mainAlignView(bamFile,rngs,rngs,filtSings=filtsings,findSplitReads=findSplitReads,dedup=dedup,filterbyMM=TRUE,MM=allowedMM[ii],returnScoreOnly=scoreOnly,
        indelRate=indelRate,mmRate=mmRate,readLength,gapOpeningArg = gapOpeningArg[ii], gapExtensionArg = gapExtensionArg[ii],substitutionMat=substitutionMat[[ii]])
        bamnames=aligned1[[1]]

        if(verbose) print(paste('primary alignment for event',ii,'done'))

        aligned2=mainAlignView(bamFile,rngs,secondrngs,filtSings=filtsings,dedup=dedup,filterbyMM=FALSE,filterbyname=TRUE,
        returnScoreOnly=scoreOnly,filternames=bamnames,indelRate=indelRate,mmRate=mmRate,readLength,gapOpeningArg = gapOpeningArg[ii],
        gapExtensionArg = gapExtensionArg[ii],substitutionMat=substitutionMat[[ii]])
        if(verbose) print(paste('secondary alignment (1 of 2) for event',ii,'done'))

        aligned3=mainAlignView(bamFile,rngs,thirdrngs,filtSings=filtsings,dedup=dedup,filterbyMM=FALSE,filterbyname=TRUE,
        returnScoreOnly=scoreOnly,filternames=bamnames,indelRate=indelRate,mmRate=mmRate,readLength,gapOpeningArg = gapOpeningArg[ii],
        gapExtensionArg = gapExtensionArg[ii],substitutionMat=substitutionMat[[ii]])
        if(verbose) print(paste('secondary alignment (2 of 2) for event',ii,'done'))

>>>>>>> fe8fe54aa12fa5f53a9270a3804121d5b15455e9
        aligned123=list(aligned1,aligned2,aligned3)

        return(aligned123)

    }

    return(alignedall)
}
