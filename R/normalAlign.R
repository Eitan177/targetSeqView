normalAlign <-
function(bamFile,typeArg ="global-local",
         substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],
         gapOpeningArg = -4, gapExtensionArg = -1,build='hg19'){
    do.call(library, list(paste("BSgenome.Hsapiens.UCSC.", build, sep='')))
    pkgname <- paste("BSgenome.Hsapiens.UCSC.", build, sep="")
    Hsapiens <- get("Hsapiens")

    bam <- getRandomBamRecords(bamFile)
    bamnames <- getBamNames(bam)
    bamseqs <- getBamSeqs(bam)
    bampos <- getBamPos(bam)
    bamchr <- getBamChr(bam)

    ## was -10 and widh=150
    refRanges <- GRanges(seqnames=bamchr,IRanges(start=bampos-200,width=500))
    refRangesList <- GenomicRanges::split(refRanges,seqnames(refRanges))
    indToRemove <- which(is.na(match(as.numeric(as.character(seqnames(refRanges))),1:22)))

    if(length(indToRemove)>0) bamseqs <- bamseqs[-indToRemove]

    prerefs <- suppressWarnings(IRanges::lapply(refRangesList[match(1:22,names(refRangesList))],
    function(x){Views(Hsapiens[[paste('chr',gsub('chr','',seqnames(x)[1]@values),sep='')]],ranges(x))}))

    refs <- do.call(c,unname(lapply(prerefs,DNAStringSet)))

    alignFnew <- pairwiseAlignment(bamseqs,refs,type=typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)
    return(alignFnew)
}
