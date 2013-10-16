alignQuick <-
function(events,dedup=TRUE,initialExpansion=0,refexpansion=400,
         indelRate,mmRate,readLength,pairlimit=2e3,gapOpeningArg = -4, gapExtensionArg = -1,
         substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],build='hg19',
         bsbuildprefix="BSgenome.Hsapiens.UCSC.",
         rngsAlign=GRanges(),conservativeContigAlign=FALSE,validChr=c(1:22,'X','Y','M'),verbose=FALSE){

    do.call(library, list(paste(bsbuildprefix, build, sep='')))
    ##pkgname <- paste("BSgenome.Hsapiens.UCSC.", build, sep="")
    genomeName <-gsub("\\.UCSC.+$","",gsub("BSgenome\\.","",bsbuildprefix))
    ii <- NULL
    if(length(gapExtensionArg)==1) gapExtensionArg=rep(gapExtensionArg,nrow(events))
    if(length(gapOpeningArg)==1) gapOpeningArg=rep(gapOpeningArg,nrow(events))
    if(length(rngsAlign)==1 & is.list(rngsAlign)) rngsAlign=rep(list(unlist(rngsAlign)),nrow(events))
    if(is.matrix(substitutionMat)){
        ss=list()
        for(ii in 1:nrow(events)) ss[[ii]]=substitutionMat
        substitutionMat=ss
    }


    alignedall=foreach(ii=1:nrow(events),.combine='list',.multicombine=TRUE) %dopar% {


        bamFile=events[ii,'Sample']

        header<-unlist(scanBamHeader(bamFile))
        if(length(grep(paste('SN',events[ii,'Chr1'],sep=':'),header))>0) chrsub <- ''
        else chrsub <- 'chr'
            rngs=GRanges(seqnames=gsub(chrsub,'',c(events[ii,'Chr1'],events[ii,'Chr2'])),
            ranges=IRanges(start=as.numeric(c(events[ii,'Start1'],events[ii,'Start2'])),
            end=as.numeric(c(events[ii,'End1'],events[ii,'End2']))))

        rngs=resize(rngs,width=1e3,fix='center')

        if(verbose) cat(paste(ifelse(ii==1,paste('Of',nrow(events),'events now working on',''),''),ii))
        if(any(is.na(match(gsub('chr','',seqnames(rngs)),validChr)))){
            print(paste('invalid chr input for event',ii,'Score will be 0'))
            aligned<-matrix(rep(1,readLength*2),ncol=2,dimnames=list(NULL,c('pvalSV','pvalContig')))
        }else{
        aligned=mainAlignQuick(bamFile,rngs,dedup=dedup,
        indelRate=indelRate,mmRate=mmRate,typeArg="global-local",readLength=readLength,pairlimit=pairlimit,gapOpeningArg = gapOpeningArg[ii], gapExtensionArg = gapExtensionArg[ii],
        substitutionMat=substitutionMat[[ii]],conservativeContigAlign=conservativeContigAlign,genomeName=genomeName,verbose=verbose)
    }

        return(aligned)

    }

    return(alignedall)
}
