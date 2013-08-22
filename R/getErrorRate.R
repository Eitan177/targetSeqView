getErrorRate <-
function(bamFile,typeArg ="global-local",recnum=1e5,
         substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],
         gapOpeningArg = -4, gapExtensionArg = -1,bsbuildprefix="BSgenome.Hsapiens.UCSC.",build='hg19'){

        do.call(library, list(paste(bsbuildprefix, build, sep='')))
        genomeName <-gsub("\\.UCSC.+$","",gsub("BSgenome\\.","",bsbuildprefix))

    ##pkgname <- paste("BSgenome.Hsapiens.UCSC.", build, sep="")
    genome <- get(genomeName)

    bam <- getRandomBamRecords(bamFile,recnum=recnum)
    bamnames <- getBamNames(bam)
    bamseqs <- getBamSeqs(bam)
    bampos <- getBamPos(bam)
    bamchr <- getBamChr(bam)
    bamstrand <- getBamStrand(bam)

    ## was -10 and widh=150
    refRanges <- GRanges(seqnames=bamchr,IRanges(start=bampos-200,width=500))
    refRangesList <- GenomicRanges::split(refRanges,seqnames(refRanges))
    indToRemove <- which(is.na(match(as.numeric(as.character(seqnames(refRanges))),1:22)))

    if(length(indToRemove)>0){
        bamseqs <- bamseqs[-indToRemove]
        bamstrand <- bamstrand[-indToRemove]
    }

    prerefs <- suppressWarnings(IRanges::lapply(refRangesList,
    function(x){Views(genome[[paste('chr',gsub('chr','',seqnames(x)[1]@values),sep='')]],ranges(x))}))

    refs <- do.call(c,unname(lapply(prerefs,DNAStringSet)))
    refs[bamstrand=='-']<-reverseComplement(refs[bamstrand=='-'])

    alignFRnew <- pairwiseAlignment(bamseqs,refs,type=typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)
        mmSummary<-mismatchSummary(pattern(alignFRnew))$position
        IndelRate<-rep(0,max(mmSummary$Position))
        tt<-table(unlist(apply(cbind(start(unlist(insertion(alignFRnew))),end((unlist(insertion(alignFRnew))))),1,function(x){x[1]:x[2]})))
        tt<-tt[as.numeric(names(tt)) %in% mmSummary$Position]
        IndelRate[as.numeric(names(tt))] <- tt/length(alignFRnew)
        tt2<-table(unlist(apply(cbind(start(unlist(deletion(alignFRnew))),end((unlist(deletion(alignFRnew))))),1,function(x){x[1]:x[2]})))##table(unlist(start(deletion(alignFRnew))))
        tt2<-table(unlist(start(deletion(alignFRnew))))##tt2[as.numeric(names(tt2)) %in% mmSummary$Position]
        IndelRate[as.numeric(names(tt2))] <- IndelRate[as.numeric(names(tt2))]+tt2/length(alignFRnew)

    return(list(mmRate=mmSummary$Probability,indelRate=IndelRate))
    }
