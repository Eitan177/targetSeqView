getBamClipped <-
    function(bamFile,rngs,readLength=100,cutoff=25){
        para <- ScanBamParam(which=rngs,what=c("mrnm","mpos","seq","pos","strand","qname","cigar"),flag=scanBamFlag(isUnmappedQuery = FALSE),
                             tag=c("MD","NM"),reverseComplement=FALSE)

        bam <- scanBam(bamFile, param=para)

        b1ind<-grep(start(rngs)[1],names(bam))
        b2ind<-(1:2)[-b1ind]
        bam1 <- bam[[b1ind]]
        bam2 <- bam[[b2ind]]
        cigar1 <- bam1$cigar
        cigar2 <- bam2$cigar
        removeMs1<-gsub('S.+','',cigar1)
        softclipped1<-rep(0,length(removeMs1))
        softclipped1[which(!grepl('[A-Z]',removeMs1))] <- as.numeric(removeMs1[which(!grepl('[A-Z]',removeMs1))])
        removeMsEnd1<-gsub('.+[MID]','',cigar1)
        softclippedEnd1<-rep(0,length(removeMsEnd1))
        softclippedEnd1[which(grepl('S',removeMsEnd1))] <- as.numeric(gsub('S','',removeMsEnd1[which(grepl('S',removeMsEnd1))]))


        removeMs2<-gsub('S.+','',cigar2)
        softclipped2<-rep(0,length(removeMs2))
        softclipped2[which(!grepl('[A-Z]',removeMs2))] <- as.numeric(removeMs2[which(!grepl('[A-Z]',removeMs2))])
        removeMsEnd2<-gsub('.+[MID]','',cigar2)
        softclippedEnd2<-rep(0,length(removeMsEnd2))
        softclippedEnd2[which(grepl('S',removeMsEnd2))] <- as.numeric(gsub('S','',removeMsEnd2[which(grepl('S',removeMsEnd2))]))

        keep<-which(softclippedEnd1>=cutoff | softclipped1 >=cutoff)
        bam1 <- lapply(bam1,function(x) {x[keep]})
        softclippedEnd1<-softclippedEnd1[keep]
        softclipped1<-softclipped1[keep]
        keep2<-which(softclippedEnd2>=cutoff | softclipped2 >=cutoff)
        bam2 <- lapply(bam2,function(x) {x[keep2]})
        softclippedEnd2<-softclippedEnd2[keep2]
        softclipped2<-softclipped2[keep2]
        inds<-list()
        for(i in 1:2){
        softclipped <- get(c('softclipped1','softclipped2')[i])
        softclippedEnd <- get(c('softclippedEnd1','softclippedEnd2')[i])
        inds[[i]]<-apply(cbind(softclipped,softclippedEnd),1,function(x){
            start=ifelse(x[1]>x[2],1,readLength-x[2]+1);
            end=ifelse(x[1]>x[2],x[1],readLength);c(start,end)})
    }

        return(list(bam1=bam1,bam2=bam2,inds1=inds[[1]],inds2=inds[[2]]))
    }

doQuickSplits <-
    function(bamFileSplits,rngsAlign,dedup=TRUE,typeArg ="global-local",substitutionMat=nucleotideSubstitutionMatrix(match = 5, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],
             pairlimit=2e3,gapOpeningArg = -100, gapExtensionArg = 0,indelRate,mmRate,readLength,conservativeContigAlign=FALSE,genomeName='Hsapiens',
             verbose=FALSE){

        bam <- getBamClipped(bamFileSplits,rngsAlign,readLength)
        if(length(bam)==0) return(list())

#        posInd<-mapply(function(x,y){which.max(c(length(x),length(y)))},pos1m[[1]],pos2m[[1]])

##        pos1mUse<-lapply(pos1m,function(x){x[posInd==1]})
##        pos2mUse<-lapply(pos2m,function(x){x[posInd==2]})

##        allmmSV<-as.vector(do.call(c,c(pos1mUse[[1]],pos2mUse[[1]])))
##        allindelSV<-c(pos1mUse[[2]],pos2mUse[[2]])

        strand1 <- as.character(bam$bam1$strand)
        strand2 <- as.character(bam$bam2$strand)
        if(length(bam$bam1$seq)>0){
        bamseqs1 <- DNAStringSet(substring(bam$bam1$seq,bam$inds1[1,],bam$inds1[2,]))
        sc1 <- IRanges(bam$inds1[1,],bam$inds1[2,])
        sc1[strand1=='-'] <- reflect(IRanges(bam$inds1[1,strand1=='-'],bam$inds1[2,strand1=='-']),IRanges(1,readLength))
        bamseqs1[strand1=='-'] <- reverseComplement(bamseqs1[strand1=='-'])
    }else{
        bamseqs1 <-bam$bam1$seq
        sc1<-IRanges()
    }
        if(length(bam$bam2$seq)>0){
        bamseqs2 <- DNAStringSet(substring(bam$bam2$seq,bam$inds2[1,],bam$inds2[2,]))
        sc2 <- IRanges(bam$inds2[1,],bam$inds2[2,])
        sc2[strand2=='-'] <- reflect(IRanges(bam$inds2[1,strand2=='-'],bam$inds2[2,strand2=='-']),IRanges(1,readLength))
        bamseqs2[strand2=='-'] <- reverseComplement(bamseqs2[strand2=='-'])
    }else{
        bamseqs2<-bam$bam2$seq
        sc2<-IRanges()
    }

        startpos1<- start(sc1)
        startpos2<- start(sc2)

### remove duplicates ###
        if(dedup){
            dups1 <- which(Biostrings::duplicated(bamseqs1))
            dups2 <- which(Biostrings::duplicated(bamseqs2))
            if(length(dups1)>0){
                    bamseqs1 <- bamseqs1[-dups1]
                    strand1 <- strand1[-dups1]
                    startpos1 <- startpos1[-dups1]
                }
            if(length(dups2)>0){
                    bamseqs2 <- bamseqs2[-dups2]
                    strand2 <- strand2[-dups2]
                    startpos2 <- startpos2[-dups2]
                }
                    #posInd <- posInd[-dups]
        }



            ## do not use all reads if coverage is very deep

            if(length(bamseqs1) > pairlimit){
                bamseqs1 <- bamseqs1[1:pairlimit]
                strand1 <- strand1[1:pairlimit]
                startpos1<-startpos1[1:pairlimit]
            }
        if(length(bamseqs2) > pairlimit){
                bamseqs2 <- bamseqs2[1:pairlimit]
                strand2 <- strand2[1:pairlimit]
                startpos2 <-startpos2[1:pairlimit]
                #posInd <- posInd[1:pairlimit]

            }


        genome <- get(genomeName)
        ref1 <- as.character(suppressWarnings(Views(genome[[paste('chr',gsub('chr','',seqnames(rngsAlign)[1]),sep='')]],ranges(rngsAlign)[1])))
        ref2 <- as.character(suppressWarnings(Views(genome[[paste('chr',gsub('chr','',seqnames(rngsAlign)[2]),sep='')]],ranges(rngsAlign)[2])))

#            bamseqs1R<-reverseComplement(bamseqs1)
#            bamseqs2R<-reverseComplement(bamseqs2)

        aligned1 <- pairwiseAlignment(bamseqs1,ref2,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)
        aligned2 <- pairwiseAlignment(bamseqs2,ref1,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)


        aligned1R <- pairwiseAlignment(reverseComplement(bamseqs1),ref2,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)
        aligned2R <- pairwiseAlignment(reverseComplement(bamseqs2),ref1,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)

        AlignBest1 <- apply(cbind(score(aligned1),score(aligned1R)),1, which.max)
        AlignBest2 <- apply(cbind(score(aligned2),score(aligned2R)),1,which.max)


        alignF1 <- aligned1[which(AlignBest1==1)]
        alignR1 <- aligned1R[which(AlignBest1==2)]
        startposF1 <- startpos1[which(AlignBest1==1)]
        startposR1 <- startpos1[which(AlignBest1==2)]
        strandF1 <- strand1[which(AlignBest1==1)]
        strandR1 <- strand1[which(AlignBest1==2)]
        alignF2 <- aligned2[which(AlignBest2==1)]
        alignR2 <- aligned2R[which(AlignBest2==2)]
        startposF2 <- startpos2[which(AlignBest2==1)]
        startposR2 <- startpos2[which(AlignBest2==2)]
        strandF2 <- strand2[which(AlignBest2==1)]
        strandR2 <- strand2[which(AlignBest2==2)]

        startposF1<-startposF1[score(alignF1)>=100]
        strandF1 <- strandF1[score(alignF1)>=100]
        alignF1<-alignF1[score(alignF1)>=100]
        startposR1<-startposR1[score(alignR1)>=100]
        strandR1 <- strandR1[score(alignR1)>=100]
        alignR1<-alignR1[score(alignR1)>=100]
        startposF2<-startposF2[score(alignF2)>=100]
        strandF2<-strandF2[score(alignF2)>=100]
        alignF2<-alignF2[score(alignF2)>=100]
        startposR2<-startposR2[score(alignR2)>=100]
        strandR2<-strandR2[score(alignR2)>=100]
        alignR2<-alignR2[score(alignR2)>=100]

        starts1tab<-table(c(start(subject(alignF1)),start(subject(alignR1))))

        starts2tab<-table(c(start(subject(alignF2)),start(subject(alignR2))))

        mults1<-as.numeric(names(starts1tab)[starts1tab>1])
        mults2<-as.numeric(names(starts2tab)[starts2tab>1])

        startposF1<-startposF1[start(subject(alignF1)) %in% mults1]
        strandF1 <- strandF1[start(subject(alignF1)) %in% mults1]
        alignF1<-alignF1[start(subject(alignF1)) %in% mults1]
        startposR1<-startposR1[start(subject(alignR1)) %in% mults1]
        strandR1 <- strandR1[start(subject(alignR1)) %in% mults1]
        alignR1<-alignR1[start(subject(alignR1)) %in% mults1]
        startposF2<-startposF2[start(subject(alignF2)) %in% mults2]
        strandF2<-strandF2[start(subject(alignF2)) %in% mults2]
        alignF2<-alignF2[start(subject(alignF2)) %in% mults2]
        startposR2<-startposR2[start(subject(alignR2)) %in% mults2]
        strandR2<-strandR2[start(subject(alignR2)) %in% mults2]
        alignR2<-alignR2[start(subject(alignR2)) %in% mults2]


        bamseqs1fornorm<-c(DNAStringSet(gsub('[^ACGT]','',pattern(alignF1))),DNAStringSet(gsub('[^ACGT]','',pattern(alignR1))))
        strand1fornorm<-c(strandF1,strandR1)
        startpos1fornorm<-c(startposF1,startposR1)
        bamseqs2fornorm<-c(DNAStringSet(gsub('[^ACGT]','',pattern(alignF2))),DNAStringSet(gsub('[^ACGT]','',pattern(alignR2))))
        strand2fornorm<-c(strandF2,strandR2)
        startpos2fornorm<-c(startposF2,startposR2)

        if(conservativeContigAlign){

            bamseqs1Rnorm<-reverseComplement(bamseqs1fornorm)
            bamseqs2Rnorm<-reverseComplement(bamseqs2fornorm)

        }else{
        bamseqs1fornorm[strand1fornorm=='-'] <- reverseComplement(bamseqs1fornorm[strand1fornorm=='-'])
        bamseqs2fornorm[strand2fornorm=='-'] <- reverseComplement(bamseqs2fornorm[strand2fornorm=='-'])
    }

        aligned1norm <- pairwiseAlignment(bamseqs1fornorm,ref1,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)
        aligned2norm <- pairwiseAlignment(bamseqs2fornorm,ref2,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)

        if(conservativeContigAlign){
        aligned1Rnorm <- pairwiseAlignment(bamseqs1Rnorm,ref1,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)
        aligned2Rnorm <- pairwiseAlignment(bamseqs2Rnorm,ref2,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)

        AlignBest1norm <- apply(cbind(score(aligned1norm),score(aligned1Rnorm)),1,which.max)
        AlignBest2norm <-apply(cbind(score(aligned2norm),score(aligned2Rnorm)),1,which.max)

        alignedF1norm <- aligned1norm[which(AlignBest1norm==1)]
        alignedR1norm <- aligned1Rnorm[which(AlignBest1norm==2)]
        alignedF2norm <- aligned2norm[which(AlignBest2norm==1)]
        alignedR2norm <- aligned2Rnorm[which(AlignBest2norm==2)]
    }else{
        alignedF1norm <- aligned1norm[strand1fornorm =='+']
        alignedR1norm <- aligned1norm[strand1fornorm =='-']
        alignedF2norm <- aligned2norm[strand2fornorm =='+']
        alignedR2norm <- aligned2norm[strand2fornorm =='-']
    }

        return(list(alignF1=alignF1,alignR1=alignR1,alignF2=alignF2,alignR2=alignR2,
                    startposF1=startposF1,startposR1=startposR1,startposF2=startposF2,startposR2=startposR2,
                    alignedF1norm=alignedF1norm,alignedR1norm=alignedR1norm,alignedF2norm=alignedF2norm,alignedR2norm=alignedR2norm))

    }
