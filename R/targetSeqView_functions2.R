##########################################################################################################
##########################################################################################################
############## FUNCTIONS TO PERFORM PARTIAL SMITH-WATERMAN ALIGNMENT AND SCORING ###########################
##########################################################################################################
##########################################################################################################


parseCigarAndMD <- function (md,cigar){

md[is.na(md)]="0"
## removed plus from [ACGT]+
mds<- gsub('[ ]$','',gsub('[ ]+',' ',gsub("([\\^]*[ACGT])[0]*", " \\1 ", md)))

removeMs<-gsub('S.+','',cigar)
softclipped<-rep(0,length(removeMs))
softclipped[which(!grepl('[A-Z]',removeMs))] <- as.numeric(removeMs[which(!grepl('[A-Z]',removeMs))])


removeMsEnd<-gsub('.+[MID]','',cigar)
softclippedEnd<-rep(0,length(removeMsEnd))
softclippedEnd[which(grepl('S',removeMsEnd))] <- as.numeric(gsub('S','',removeMsEnd[which(grepl('S',removeMsEnd))]))



firstM <- as.numeric(gsub('.+[A-Z]','',gsub('M.*$','',cigar)))
if(length(grep('[ID]',cigar))==0) {
    bonifiedInsOrDels <- rep(0,length(cigar))
}else{
bonifiedInsOrDels<-unlist(lapply(vmatchPattern('D',cigar),length))+
        unlist(lapply(vmatchPattern('I',cigar),length))
}


indels<-apply(cbind(softclipped,bonifiedInsOrDels,firstM,softclippedEnd),1,function(x){list(c(ifelse(x[1]>0,1,0),
                                                                                 ifelse(x[2] > 0,(x[1]+x[3]),0),
                                                                                 ifelse(x[4]>0,x[4],0)))})

#indelsInd<-regexpr('[IDS]',cigar)
#indels<-rep(0,length(cigar))
#indels[indelsInd != -1] = indelsInd[indelsInd != -1]


toag<-rep(seq_along(mds),vcountPattern(' ',mds)+1)
pp<-paste(mds,collapse=' ')
ss<-strsplit(pp,'[ ]+')[[1]]
ss[grep('\\^',ss)]=nchar(grep('\\^',ss,value=TRUE))-1
ss[ss %in% c('A','C','G','T')]=.99

    tosee <- aggregate(ss,by=list(toag),function(x){cs=cumsum(as.numeric(as.character(x)));cs[x==.99]},simplify=FALSE)

if(is.matrix(tosee[,2])) agpos <- do.call(c,apply(tosee[,2],1,list))
else agpos <- tosee[,2]


mmpositions<-mapply(function(x,y){x+y},agpos,softclipped,SIMPLIFY=FALSE)

return(list(mmpositions,indels,softclip=cbind(softclipped,softclippedEnd)))

}

getMismatchQuick<-function(readLength,alignFmm,alignRmm,softclippedseq){

        mmtableF<-mismatchTable(alignFmm)
        mmtableF$PatternStart[mmtableF$PatternStart > readLength] <- readLength
        softclippedF<-softclippedseq[mmtableF$PatternId,,drop=FALSE]

        ## we don't want to penalize a read for having a mismatch in our realignment in an area that was softclipped in the bam
        mmtableF<-mmtableF[mmtableF$PatternStart > softclippedF[,'softclipped'] & mmtableF$PatternStart < softclippedF[,'softclippedEnd'],]
        forwardMtable <- table(mmtableF$PatternStart)

        mmtableR<-mismatchTable(alignRmm)
        mmtableR$PatternStart[mmtableR$PatternStart > readLength] <- readLength
        mmtableR$correctedPStart<-readLength-mmtableR$PatternStart+1
        softclippedR<-softclippedseq[mmtableR$PatternId,,drop=FALSE]

        ## we don't want to penalize a read for having a mismatch in our realignment in an area that was softclipped in the bam
        mmtableR<-mmtableR[mmtableR$PatternStart > softclippedR[,'softclipped'] & mmtableR$PatternStart < softclippedR[,'softclippedEnd'],]
        reverseMtable <- table(mmtableR$correctedPStart)

        allMismatches <- rep(0,readLength)

        allMismatches[as.numeric(names(forwardMtable))] <- forwardMtable
        allMismatches[as.numeric(names(reverseMtable))]<- allMismatches[as.numeric(names(reverseMtable))]+reverseMtable

        return(allMismatches)
    }


getindelsQuick<-function(readLength,alignFmm,alignRmm,softclippedseq){

        allIndel <- rep(0,readLength)

        insertionsF<-do.call(c,mapply(function(x,y){setdiff(reduce(x),y);reduce(x)},insertion(alignFmm),deletion(alignFmm)))
        if(is.null(insertionsF)) insertionsF <- IRanges()
        width(insertionsF)[width(insertionsF)>25]=25

        insertionsR<-do.call(c,mapply(function(x,y){setdiff(reduce(x),y);reduce(x)},insertion(alignRmm),deletion(alignRmm)))

        if(is.null(insertionsR)) insertionsR <- IRanges()
        width(insertionsR)[width(insertionsR)>25]=25
        insertionInds<-unlist(apply(cbind(start(insertionsF),(end(insertionsF))),1,function(x){x[1]:x[2]}))
        ### just for deletion positions

        deletionInds<-unlist(apply(cbind(start(unlist(deletion(alignFmm))),end(unlist(deletion(alignFmm)))),1,function(x){x[1]:x[2]}))

        delf<-unlist(deletion(alignFmm))
##        delf<-delf[refInd[start(delf)] == refInd[end(delf)]]

        forwardIndeltable<-table(c(insertionInds,start(delf)))
### just for deletions positions
        forwardIndeltablewDel<-table(c(insertionInds,deletionInds))
        forwardIndeltablewDel<-forwardIndeltablewDel[as.numeric(names(forwardIndeltablewDel)) %in% 1:readLength]

        forwardIndeltable<-forwardIndeltable[as.numeric(names(forwardIndeltable)) %in% 1:readLength]

        insertionIndsR<-unlist(apply(cbind(start(insertionsR),(end(insertionsR))),1,function(x){x[1]:x[2]}))
## just for deletions positions
        deletionIndsR<-unlist(apply(cbind(start(unlist(deletion(alignRmm))),(end(unlist(deletion(alignRmm))))),1,function(x){x[1]:x[2]}))
        reverseIndeltable<-table(readLength-c(insertionIndsR,start(unlist(deletion(alignRmm))))+1)
## just for deletions positions
        reverseIndeltableWdel<-table(readLength-c(insertionIndsR,deletionIndsR)+1)
        reverseIndeltableWdel<-reverseIndeltableWdel[as.numeric(names(reverseIndeltableWdel)) %in% 1:readLength]

        reverseIndeltable<-reverseIndeltable[as.numeric(names(reverseIndeltable)) %in% 1:readLength]
        allIndel[as.numeric(names(forwardIndeltable))]<- forwardIndeltable
        allIndel[as.numeric(names(reverseIndeltable))]<- allIndel[as.numeric(names(reverseIndeltable))]+reverseIndeltable

allindelsWdel<-rep(0,readLength)
allindelsWdel[as.numeric(names(forwardIndeltablewDel))]<- forwardIndeltablewDel
allindelsWdel[as.numeric(names(reverseIndeltableWdel))]<- allindelsWdel[as.numeric(names(reverseIndeltableWdel))]+reverseIndeltableWdel

        return(list(allIndel,allindelsWdel))
    }

getRandomBamRecords <-
    function(bamFile,reverseC=TRUE,recnum=1e5){
        bamheader<-scanBamHeader(bamFile)[[1]]
        chr<- names(bamheader$targets)
        chrlen <-bamheader$targets
        rngsCheck<-unlist(GRangesList(sapply(1:length(chr),function(x){
            GRanges(seqnames=chr[x],ranges= breakInChunks(chrlen[x],100e3))})))
            ind <- nrecs<-0
            para <- ScanBamParam()
        while(nrecs < recnum){
            ind <- ind+1
        para <- ScanBamParam(what=scanBamWhat(),reverseComplement=reverseC,
                             which=rngsCheck[1:ind],
                             flag=scanBamFlag(isProperPair=TRUE))
        cb<-countBam(bamFile,param=para)

            nrecs<-sum(cb$records)
        print(paste('found',nrecs,'normal reads, looking for another', max(c(0,recnum-nrecs))))
        }
        bam <- scanBam(bamFile,param=para)
        return(bam)
    }

getBam <-
    function(bamFile,rngs){
        para <- ScanBamParam(which=rngs,what=c("mrnm","mpos","seq","pos","strand","qname","cigar"),flag=scanBamFlag(isUnmappedQuery = FALSE),
                             tag=c("MD","NM"),reverseComplement=TRUE)

        bam <- scanBam(bamFile, param=para)


        b1ind<-grep(start(rngs)[1],names(bam))
        b2ind<-(1:2)[-b1ind]
        bam1 <- bam[[b1ind]]
        bam1$MD=bam1$tag$MD
        bam2 <- bam[[b2ind]]
        bam2$MD=bam2$tag$MD

        keep<-which(!(is.na(bam1$mpos)))
        bam1 <- lapply(bam1,function(x) {x[keep]})

        if(length(bam1$seq)==0) return (list())

        mateRanges1<- GRanges(seqnames=as.character(bam1$mrnm),IRanges(start=as.numeric(bam1$mpos),width=1))
        ovskeep<-which(!is.na(findOverlaps(mateRanges1,rngs[2],select="first")))#,maxgap=300)))

        ## if we don't have any good mates, keep a few anyway to test
        if(length(ovskeep)==0) {
            qnametab <-table(bam1$qname)
            qnamekeep<-names(qnametab[qnametab ==2])
            if(length(qnamekeep)==0) return(list())
            bam1ind<-which(bam1$qname %in% c(qnamekeep))
            sqnamebam1<-split(bam1ind,bam1$qname[bam1ind])
            b1ind<-unname(unlist(lapply(sqnamebam1,function(x) {x[1]})))
            b2ind<-unname(unlist(lapply(sqnamebam1,function(x) {x[2]})))
            bam1filt <-lapply(bam1,function(x){x[b1ind]})
            bam2filt <-lapply(bam1,function(x){x[b2ind]})
        }else{
        bam1filt<-lapply(bam1,function(x){x[ovskeep]})
        keep2<-match(bam1filt$qname,bam2$qname)
        bam2filt<-lapply(bam2,function(x){x[keep2]})
    }

        return(list(bam1filt,bam2filt))
    }

precomputedTargetCapture100bpMMRate<-
    function() {
        return(c(0.0114,0.0084,0.0096,0.0083,0.0076,0.0069,0.006,0.0063,0.0049,0.0047,0.004,0.0042,0.0035,0.0034,0.0034,0.0034,0.0031,0.0033,0.0036,0.0033,0.0031,0.0032,0.0034,0.0036,0.0032,0.0038,0.0033,0.0037,0.0036,0.0038,0.0034,0.0034,0.0034,0.0036,0.0038,0.0038,0.0037,0.0036,0.0032,0.0039,0.0037,0.0035,0.0035,0.0035,0.0041,0.0044,0.0042,0.0039,0.0042,0.0045,0.0038,0.0041,0.0043,0.0043,0.0049,0.0048,0.005,0.0048,0.0054,0.0052,0.0059,0.0062,0.0066,0.0065,0.0067,0.0074,0.0075,0.0078,0.0084,0.0079,0.0081,0.008,0.0085,0.0095,0.0096,0.0092,0.0094,0.01,0.011,0.0112,0.0118,0.0122,0.0135,0.0128,0.0143,0.0149,0.0161,0.0164,0.0168,0.018,0.0185,0.02,0.021,0.0218,0.0219,0.024,0.0245,0.0252,0.0284,0.0525))
    }
precomputedTargetCapture100bpIndelRate<-
    function(){
        return(c(0.01099,0.01012,0.0061,0.0061,0.0054,0.00451,0.00353,0.00297,0.00253,0.00209,0.00153,0.00136,0.00117,0.00096,0.00079,0.00097,8e-04,0.00062,7e-04,0.00065,0.00063,7e-04,0.00049,5e-04,0.00054,5e-04,0.00056,5e-04,0.00045,0.00049,0.00025,0.00045,0.00035,5e-04,0.00033,0.00041,0.00053,4e-04,0.00044,0.00031,0.00037,0.00045,0.00043,0.00047,0.00034,3e-04,0.00035,0.00032,0.00027,0.00033,0.00045,0.00038,0.00044,0.00043,0.00035,0.00042,0.00038,0.00043,4e-04,0.00036,0.00037,0.00038,0.00033,0.00034,0.00064,0.00042,5e-04,4e-04,0.00039,0.00042,0.00047,0.00086,0.00045,0.00039,5e-04,0.00046,0.00041,0.00047,0.00044,0.00055,0.00071,0.00062,0.00062,0.00073,0.00074,0.00058,0.00073,0.00078,0.00067,0.00087,0.00103,0.00093,0.00126,0.00135,0.00129,0.00166,0.00195,0.00197,0.00226,6e-05))
    }

precomputedWholeGenome101bpMMRate <-
    function(){
        return(c(0.0056,0.0065,0.0048,0.0047,0.0056,0.0056,0.0061,0.0059,0.0053,0.0052,0.0064,0.0059,0.005,0.0061,0.0063,0.0061,0.0065,0.006,0.0062,0.0059,0.0062,0.0059,0.0062,0.0065,0.0064,0.0062,0.0066,0.0066,0.0066,0.0074,0.0074,0.0073,0.0071,0.0082,0.0099,0.0099,0.0088,0.0097,0.0107,0.0129,0.0113,0.0107,0.011,0.0116,0.0116,0.0125,0.0156,0.0135,0.0133,0.0126,0.0112,0.0122,0.0116,0.0118,0.0143,0.0117,0.0125,0.012,0.0131,0.013,0.0144,0.0133,0.0152,0.0137,0.0148,0.0157,0.0142,0.0149,0.0161,0.0169,0.0155,0.0172,0.0168,0.0174,0.0171,0.0158,0.0171,0.0175,0.0168,0.0163,0.0178,0.0178,0.0174,0.0184,0.0191,0.0198,0.0197,0.0199,0.0203,0.0195,0.0204,0.0201,0.0206,0.0209,0.0249,0.0227,0.0224,0.0227,0.0242,0.0224,0.0273))
    }
precomputedWholeGenome101bpIndelRate <-
    function(){
        return(c(0.00022,0.00039,0.00046,0.00039,0.00043,0.00044,0.00052,0.00046,0.00036,0.00049,5e-04,0.00043,0.00045,0.00053,0.00048,0.00036,0.00043,0.00032,0.00036,0.00039,0.00038,0.00042,0.00038,0.00038,0.00038,0.00048,0.00041,0.00046,0.00038,0.00031,0.00041,0.00036,0.00044,4e-04,6e-04,0.00063,7e-04,0.00075,0.00087,0.00078,0.00087,0.00074,0.00084,0.00095,0.00093,0.00091,0.00101,0.00095,0.00097,0.00095,0.00093,9e-04,0.00081,0.00097,0.00104,0.0011,0.0011,0.00111,0.00114,0.00129,0.00131,0.00133,0.00135,0.00144,0.00159,0.00158,0.00153,0.00145,0.00138,0.00158,0.00177,0.00165,0.00186,0.00192,0.00162,0.00196,0.00172,0.00188,0.00191,0.00203,0.00188,0.00199,0.00205,0.00209,0.00201,0.00222,0.00214,0.00247,0.00258,0.0026,0.00278,0.00289,0.0027,0.00281,0.0028,0.00289,0.00284,0.00338,0.00306,0.00252,5e-04))
    }

mainAlignQuick <-    function(bamFile,rngsAlign,dedup=TRUE,typeArg ="global-local",substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],pairlimit=2e3,gapOpeningArg = -4, gapExtensionArg = -1,indelRate,mmRate,readLength=100,conservativeContigAlign=FALSE,genomeName='Hsapiens',verbose=FALSE){

     alignedformatted <- doAlignAndformatQuick(bamFile,rngsAlign,dedup,typeArg,substitutionMat,pairlimit,gapOpeningArg ,gapExtensionArg,indelRate=indelRate,mmRate=mmRate,
                                               readLength=readLength,conservativeContigAlign=conservativeContigAlign,genomeName=genomeName,verbose=verbose)

     if(length(alignedformatted)==0) return(matrix(rep(1,readLength*2),ncol=2,dimnames=list(NULL,c('pvalSV','pvalContig'))))

    alignF1<-alignedformatted[[1]]
    alignR1<-alignedformatted[[2]]
    alignF2<-alignedformatted[[3]]
    alignR2<-alignedformatted[[4]]
     allmmSV<-alignedformatted[[5]]
     allindelSV<-alignedformatted[[6]]
     softclippedseq1<-alignedformatted$sf1
     softclippedseq2<-alignedformatted$sf2

     pvalside1<-getpvalQuick(alignF1,alignR1,readLength,mmRate,indelRate,softclippedseq1)
     pvalside2<-getpvalQuick(alignF2,alignR2,readLength,mmRate,indelRate,softclippedseq2)

     pvalContig <- pvalside1*pvalside2

     pvalSV <- getpvalQuickSV(allmmSV,allindelSV,readLength,length(alignF1)+length(alignF2)+length(alignR1)+length(alignR2),mmRate,indelRate)

     return(cbind(pvalSV,pvalContig))


 }


getpvalQuickSV<-function(allmmSV,allindelSV,readLength,numreads,mmRate,indelRate){

    tt <-table(allmmSV)
    allMismatches <- rep(0,readLength)
    allMismatches[ceiling(as.numeric(names(tt)))]<-tt

        allIndel<-rep(0,readLength)
    allindeltt<-table(unname(unlist(allindelSV)))
    allindeltt <- allindeltt[names(allindeltt)!= '0']
    allIndel[as.numeric(names(allindeltt))]<-allindeltt

    #assignIndels<-table(sample(1:readLength,sum(allindelSV),replace=TRUE))
    #allIndel[as.numeric(names(assignIndels))]<-assignIndels



        mismatchBinom <- sapply(1:readLength,function(x){sum(dbinom(allMismatches[x]:numreads,numreads,mmRate[x]))})
        indelBinom <- sapply(1:readLength,function(x){sum(dbinom(allIndel[x]:numreads,numreads,indelRate[x]))})
        pval <- indelBinom*mismatchBinom ##c(indelBinom,mismatchBinom)
        return(pval)
}


getpvalQuick<-function(alignF,alignR,readLength,mmRate,indelRate,softclippedseq){

        ### new
        allMismatches<-getMismatchQuick(readLength,alignF,alignR,softclippedseq)
        allIndel<-getindelsQuick(readLength,alignF,alignR,softclippedseq)
        allIndel<-allIndel[[1]]


        numreads<-length(alignF)+length(alignR)
        mismatchBinom <- sapply(1:readLength,function(x){sum(dbinom(allMismatches[x]:numreads,numreads,mmRate[x]))})

        indelBinom <- sapply(1:readLength,function(x){sum(dbinom(allIndel[x]:numreads,numreads,indelRate[x]))})

        pval <- indelBinom*mismatchBinom ##c(indelBinom,mismatchBinom)
        return(pval)
}


doAlignAndformatQuick <-
    function(bamFile,rngsAlign,dedup=TRUE,typeArg ="global-local",substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],
             pairlimit=2e3,gapOpeningArg = -4, gapExtensionArg = -1,indelRate,mmRate,readLength,conservativeContigAlign=FALSE,genomeName='Hsapiens',verbose=FALSE){

        bam <- getBam(bamFile,rngsAlign)
        if(length(bam)==0) return(list())

        pos1m<-parseCigarAndMD(bam[[1]]$MD,bam[[1]]$cigar)
        pos2m<-parseCigarAndMD(bam[[2]]$MD,bam[[2]]$cigar)
        posInd<-mapply(function(x,y){which.max(c(length(x),length(y)))},pos1m[[1]],pos2m[[1]])

        pos1mUse<-lapply(pos1m,function(x){x[posInd==1]})
        pos2mUse<-lapply(pos2m,function(x){x[posInd==2]})

        allmmSV<-as.vector(do.call(c,c(pos1mUse[[1]],pos2mUse[[1]])))
        allindelSV<-c(pos1mUse[[2]],pos2mUse[[2]])


        bamseqs1 <- bam[[1]]$seq
        bamseqs2 <- bam[[2]]$seq
        matestrand1 <- as.character(bam[[2]]$strand)
        matestrand2 <- as.character(bam[[1]]$strand)

        ### remove duplicates ###
        if(dedup){
            dups1 <- which(Biostrings::duplicated(bamseqs1))
            dups2 <- which(Biostrings::duplicated(bamseqs1))
            if(length(dups1)>0 & length(dups2) >0){
                dupstab <- table(c(dups1,dups2))
                dupstab2<- dupstab[dupstab==2]
                if(length(dupstab2)>0){
                    dups <- as.numeric(names(dupstab2))
                    bamseqs1 <- bamseqs1[-dups]
                    bamseqs2 <- bamseqs2[-dups]
                    matestrand1 <- matestrand1[-dups]
                    matestrand2 <- matestrand2[-dups]
                    posInd <- posInd[-dups]
                }
            }
        }

            ## do not use all reads if coverage is very deep
            if(length(bamseqs1) > pairlimit){

                bamseqs1 <- bamseqs1[1:pairlimit]
                bamseqs2 <- bamseqs2[1:pairlimit]
                matestrand1 <- matestrand1[1:pairlimit]
                matestrand2 <- matestrand2[1:pairlimit]
                posInd <- posInd[1:pairlimit]

            }


        genome <- get(genomeName)
        ref1 <- as.character(suppressWarnings(Views(genome[[paste('chr',gsub('chr','',seqnames(rngsAlign)[1]),sep='')]],ranges(rngsAlign)[1])))
        ref2 <- as.character(suppressWarnings(Views(genome[[paste('chr',gsub('chr','',seqnames(rngsAlign)[2]),sep='')]],ranges(rngsAlign)[2])))

        ##revert back
        if(conservativeContigAlign){
        ## might not need this
            bamseqs1R<-reverseComplement(bamseqs1)
            bamseqs2R<-reverseComplement(bamseqs2)
        }else{
        bamseqs1[matestrand1=='+'] <- reverseComplement(bamseqs1[matestrand1=='+'])
        bamseqs2[matestrand2=='+'] <- reverseComplement(bamseqs2[matestrand2=='+'])
    }

        aligned1 <- pairwiseAlignment(bamseqs1,ref2,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)
        aligned2 <- pairwiseAlignment(bamseqs2,ref1,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)

        if(conservativeContigAlign){
        ## might not need this
        aligned1R <- pairwiseAlignment(bamseqs1R,ref2,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)
        aligned2R <- pairwiseAlignment(bamseqs2R,ref1,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)
        AlignBest <- apply(cbind(score(aligned1),score(aligned1R),score(aligned2),score(aligned2R)),1,which.max)
        # might not need this
        alignedF1 <- aligned1[which(AlignBest==1)]
        alignedR1 <- aligned1R[which(AlignBest==2)]
        alignedF2 <- aligned2[which(AlignBest==3)]
        alignedR2 <- aligned2R[which(AlignBest==4)]

    }else{
        AlignBest <- apply(cbind(score(aligned1),score(aligned2)),1,which.max)
        alignedF1 <- aligned1[which(AlignBest==1 & matestrand1 !='+')]
        alignedR1 <- aligned1[which(AlignBest==1 & matestrand1 =='+')]
        alignedF2 <- aligned2[which(AlignBest==2 & matestrand2 !='+')]
        alignedR2 <- aligned2[which(AlignBest==2 & matestrand2 =='+')]

    }

        return(list(alignedF1,alignedR1,alignedF2,alignedR2,allmmSV,allindelSV,sf1=pos1m$softclip,sf2=pos2m$softclip))

    }

##########################################################################################################
##########################################################################################################
############## FUNCTIONS TO PERFORM FULL SMITH-WATERMAN ALIGNMENT, SCORING, AND VIEWING ##################
##########################################################################################################
##########################################################################################################

getMismatchFull<-function(readLength,alignFmm,alignRmm){

        mmtableF<-mismatchTable(alignFmm)
        mmtableF$PatternStart[mmtableF$PatternStart > readLength] <- readLength
        forwardMtable <- table(mmtableF$PatternStart)

        mmtableR<-mismatchTable(alignRmm)
        mmtableR$PatternStart[mmtableR$PatternStart > readLength] <- readLength
        mmtableR$correctedPStart<-readLength-mmtableR$PatternStart+1
        reverseMtable <- table(mmtableR$correctedPStart)

        allMismatches <- rep(0,readLength)

        allMismatches[as.numeric(names(forwardMtable))] <- forwardMtable
        allMismatches[as.numeric(names(reverseMtable))]<- allMismatches[as.numeric(names(reverseMtable))]+reverseMtable
        return(allMismatches)
    }


getindelsFull<-function(readLength,alignFmm,alignRmm,refInd){

        allIndel <- rep(0,readLength)

        insertionsF<-do.call(c,mapply(function(x,y){setdiff(reduce(x),y);reduce(x)},insertion(alignFmm),deletion(alignFmm)))
        if(is.null(insertionsF)) insertionsF <- IRanges()
        width(insertionsF)[width(insertionsF)>25]=25

        insertionsR<-do.call(c,mapply(function(x,y){setdiff(reduce(x),y);reduce(x)},insertion(alignRmm),deletion(alignRmm)))

        if(is.null(insertionsR)) insertionsR <- IRanges()
        width(insertionsR)[width(insertionsR)>25]=25
        insertionInds<-unlist(apply(cbind(start(insertionsF),(end(insertionsF))),1,function(x){x[1]:x[2]}))

        deletionInds<-unlist(apply(cbind(start(unlist(deletion(alignFmm))),end(unlist(deletion(alignFmm)))),1,function(x){x[1]:x[2]}))

        delf<-unlist(deletion(alignFmm))
        delf<-delf[refInd[start(delf)] == refInd[end(delf)]]

        forwardIndeltable<-table(c(insertionInds,start(delf)))

        forwardIndeltablewDel<-table(c(insertionInds,deletionInds))
        forwardIndeltablewDel<-forwardIndeltablewDel[as.numeric(names(forwardIndeltablewDel)) %in% 1:readLength]

        forwardIndeltable<-forwardIndeltable[as.numeric(names(forwardIndeltable)) %in% 1:readLength]

        insertionIndsR<-unlist(apply(cbind(start(insertionsR),(end(insertionsR))),1,function(x){x[1]:x[2]}))
## just for deletions positions
        deletionIndsR<-unlist(apply(cbind(start(unlist(deletion(alignRmm))),(end(unlist(deletion(alignRmm))))),1,function(x){x[1]:x[2]}))
        reverseIndeltable<-table(readLength-c(insertionIndsR,start(unlist(deletion(alignRmm))))+1)
## just for deletions positions
        reverseIndeltableWdel<-table(readLength-c(insertionIndsR,deletionIndsR)+1)
        reverseIndeltableWdel<-reverseIndeltableWdel[as.numeric(names(reverseIndeltableWdel)) %in% 1:readLength]

        reverseIndeltable<-reverseIndeltable[as.numeric(names(reverseIndeltable)) %in% 1:readLength]
        allIndel[as.numeric(names(forwardIndeltable))]<- forwardIndeltable
        allIndel[as.numeric(names(reverseIndeltable))]<- allIndel[as.numeric(names(reverseIndeltable))]+reverseIndeltable

allindelsWdel<-rep(0,readLength)
allindelsWdel[as.numeric(names(forwardIndeltablewDel))]<- forwardIndeltablewDel
allindelsWdel[as.numeric(names(reverseIndeltableWdel))]<- allindelsWdel[as.numeric(names(reverseIndeltableWdel))]+reverseIndeltableWdel

        return(list(allIndel,allindelsWdel))
    }


dfmake2 <-
    function(fordf1A,deL,bamn,refalign,refInd,refRa,refM,filtSings,bamnamesRemove,findSplitReads,rngsAlign,doingSplits,pairedo){
        ## combine paired end

        fordPEs <- split(fordf1A,bamn)

        if(filtSings){
            ## make sure two reads pass per read name pass to this point
            if(findSplitReads == FALSE){
                btab <- table(bamn);fordPEs <- fordPEs[names(fordPEs) %in% names(btab)[btab>1]]
            }

            ## make sure for each pair of read, majority of reference is on different reference if there are two reference

            if(length(rngsAlign)>1 & length(fordPEs)>1){
                difrefs<-FALSE
                splitlim=30
                ## added & splitlim>0
                while(all(difrefs == FALSE) & splitlim>0){
                difrefs <- sapply(fordPEs,function(x){
                    userec=FALSE;
                    if(findSplitReads){
                        freadref=table(x[[1]][,3]); SplitRead <- length(freadref)>1 & all(freadref>splitlim)
                        if(length(x)>1){
                            sreadref=table(x[[2]][,3]); SecondReadSplit <- length(sreadref)>1 & all(sreadref>splitlim)
                            SplitRead <- SplitRead | SecondReadSplit
                        }
                        userec <- SplitRead

                    };
                    if((userec==FALSE | findSplitReads == FALSE) & length(x)>1){
                        freadref=table(x[[1]][,3]);sreadref=table(x[[2]][,3]);
                        userec <- names(freadref)[which.max(freadref)] != names(sreadref)[which.max(sreadref)]

                    }
                    userec
                })
                splitlim=splitlim-5
#                print(splitlim)
                if(splitlim<6){print('splitlim less than 6')
                           print('splitlim less than 6')
            print('splitlim less than 6')
            #print(rngsAlign)
                           }
            }

                if(!any(difrefs)){
                    hack1<-which(unlist(unname(sapply(fordPEs,function(x){length(table(do.call(rbind,x)[,3]))==2}))))[1]
                    difrefs <-ifelse(!is.na(hack1)>0,hack1,1)
                    if(doingSplits) return(data.frame())
                }
                fordPEs <- fordPEs[difrefs]
            }
        }


        fordPEsComb <- lapply(fordPEs,function(x){
            if(length(x)<=2){
                do.call(rbind,x)
            }else{
                if(length(rngsAlign)>1){
                    readrefs=unlist(lapply(x,function(y){y[1,3]}))
                    tokeep=match(unique(readrefs),readrefs)
                }else{
                    di=as.matrix(dist(unlist(lapply(x,function(y){y[1,2]}))))
                    tokeep=which(di == max(di), arr.ind = TRUE)[1:2]
                }
                do.call(rbind,x[tokeep])
            }
        })
        if(length(fordPEsComb)>1){
            ### added for findSplitReadsOnly
        fordPEsCombDup<- list();readfrac<-.7
        while(length(fordPEsCombDup)==0){
            fordPEsCombDup <- lapply(fordPEsComb,function(x){x2=x[order(as.numeric(x[,2]),x[,1]),];x3=x2[!duplicated(x2[,2]),];
                                                             if(length(which(duplicated(x2[,2])))>(nrow(x3)*readfrac) & length(rngsAlign)>1){NA}else{x3}})
            fordPEsCombDup[which(is.na(fordPEsCombDup))]<-NULL
            if(readfrac==1) fordPEsCombDup <- fordPEsComb[1]
            readfrac<-readfrac+.1
        }
        ## end added for findSplitReadsOnly
    }else fordPEsCombDup <- fordPEsComb

        if(length(rngsAlign)>1 & (bamnamesRemove[1]=='' | doingSplits==TRUE)){
        orientation<-unlist(unname(lapply(fordPEsCombDup,function(x){

            if(length(reduce(matchPattern('E',paste(x[,1],collapse=''))))==1){
                o<-'a.facing'

            }else if(length(reduce(matchPattern(c('M'),gsub('[ACGT]','M',paste(x[,1],collapse='')))))==1){
                o<-'c.back2back'
            }else{
            o<-'b.probableInv'
        }
            ### added | pairedo=='' for findSplitReadsOnly
            gaploc<-reduce(matchPattern('-',paste(x[,1],collapse='')))
            nocontradiction <- ((length(gaploc)==2 & !any(is.na(match(c('a.facing','c.back2back'),c(o,pairedo))))) |
            length(gaploc)==1 & (o==pairedo | pairedo==''))

            if(doingSplits & (length(rngsAlign)==1 | (nocontradiction &
               (length(gaploc)==1 | (length(gaploc)==2 & !any(is.na(match(c(1,nrow(x)), c(start(gaploc),end(gaploc)))))))))){
                o<-'d.split'
            }

            mistakes<-length(matchPattern(c('N'),gsub('[ACGT]','N',paste(x[,1],collapse=''))))
            ret<-mistakes;names(ret)<-o;return(ret)

        })))

        splitOrientation<-split(orientation,names(orientation))
        meanMistakes<-sapply(splitOrientation,mean)
        ##readsnum<-sapply(splitOrientation,length)
        typeuse<-names(meanMistakes)[which.min(meanMistakes)]
        ##typeuse<-names(meanMistakes)[which.max(readsnum)]

        ##typeuse<-names(table(orientation))[which.max(table(orientation))]
        ##if(typeuse != 'c.back2back') fordPEsCombDup<-fordPEsCombDup[orientation != 'c.back2back']

        if(doingSplits) {

            fordPEsCombDup<-fordPEsCombDup[names(orientation) == 'd.split']
            ### removed typeuse for findSplitReadsOnly
            if(length(fordPEsCombDup)==0) return(data.frame())
        }else if(length(which(names(orientation)==typeuse)) >2) fordPEsCombDup<-fordPEsCombDup[names(orientation) ==typeuse]

    }

        vvA <- unlist(lapply(fordPEsCombDup,function(x){as.numeric(as.character(data.frame(x)[1,2]))}))

#if(length(rngsAlign)==1 & doingSplits)browser()
        fordf2A <- sapply(order(as.numeric(vvA)),function(x){list(fordPEsCombDup[[x]])})

        fordf3A <- sapply(1:length(fordf2A),function(x){list(cbind(fordf2A[[x]],rep((x),nrow(fordf2A[[x]]))))})


        names(fordf3A) <- names(fordPEsCombDup)[(order(as.numeric(vvA)))]
        namecol <- rep(names(fordf3A),unlist(sapply(fordf3A,nrow)))
        dftailA <- data.frame(do.call(rbind,fordf3A),namecol)

        for(ii in c(2:5)){
            if(ii==4){chrcol <- gsub('chr','',as.character(dftailA[,ii]))
                      dftailA[,ii] <- chrcol
                  }else{
                      dftailA[,ii] <- as.numeric(gsub('chr','',as.character(dftailA[,ii])))
                  }
        }

        colnames(dftailA) <- c('afvec2','posvec2','ref','chr','ypostail','refname')

        prerefrange <- aggregate(dftailA$posvec2,by=list(dftailA$ref),FUN=range)
        prerefsee <- substring(refalign,match(prerefrange$x[,1],refInd),match(prerefrange$x[,2],refInd))
        indr <- eval(parse(text=paste('c(',paste(apply(cbind(match(prerefrange$x[,1],refInd),match(prerefrange$x[,2],refInd)),1,
                           function(x){paste(x,collapse=':')}),collapse=','),')')))
        refIndsee <- refInd[indr]
        refSide <- refM[indr]

        refChr<-gsub('chr','',as.character(seqnames(refRa)))[indr]
        refsee <- Biostrings::unlist(DNAStringSet(prerefsee))

        formatrefsee <- data.frame(t(Biostrings::as.matrix(DNAStringSet(refsee))),refIndsee,refSide,refChr, -1,'refseq')
        colnames(formatrefsee) <- colnames(dftailA)
        dftailA <- rbind(dftailA,formatrefsee)#,zeros)
        if(length(rngsAlign)>1) rownames(dftailA)[1]<-typeuse
        return(dftailA)
    }
dfmake <-
    function(aa,frvec,fname,amatname,startpos,deL,bamn,refalign,refInd,refRa,refM,filtSings,bamnamesRemove,findSplitReads,rngsAlign,doingSplits,pairedo){

        ees <- rep('E',7)
        fordf1A <- mapply(function(x,y,z){x2=strsplit(gsub('S','N',toString(x)),'')[[1]];if(z=='F'){newy=y; x2=c(x2,ees) }else{newy=max(1,y-7); x2=c(rep('E',y-newy),x2)};
                                          ## na.omit is to treat condition where ees run off the end
                                          cc=na.omit(cbind(x2,refInd[newy:(length(x2)+newy-1)],refM[newy:(length(x2)+newy-1)],
                                          as.character(seqnames(refRa))[newy:(length(x2)+newy-1)]));dd=diff(as.numeric(cc[,2]));
                                          ll=list(cc)},Biostrings::as.list(aa),as.list(startpos),as.list(frvec))
        return(dfmake2(fordf1A,deL,bamn,refalign,refInd,refRa,refM,filtSings,bamnamesRemove,findSplitReads,rngsAlign,doingSplits,pairedo))
    }
doPlot <-
    function(trimA){
        vplayout <- function(x, y){viewport(layout.pos.row = x, layout.pos.col = y)}

        grid.newpage()
        pushViewport(viewport(layout = grid.layout(1, 1)))
        print(trimA,vp = vplayout(1,1))
    }
getBamChr <-
    function(bam){
        bamChr <- vector()
        for(ii in 1:length(bam)){
            bamChr <- c(bamChr,bam[[ii]]$rname)
        }
        return(bamChr)
    }
getBamNames <-
    function(bam){
        bamnames <- vector()
        for(ii in 1:length(bam)){
            bamnames <- c(bamnames,bam[[ii]]$qname)
        }
        return(bamnames)
    }
getBamPos <-
    function(bam){
        bamPos <- vector()
        for(ii in 1:length(bam)){
            bamPos <- c(bamPos,bam[[ii]]$pos)
        }
        return(bamPos)
    }
getBamStrand <-
    function(bam){
        bamStrand <- vector()
        for(ii in 1:length(bam)){
            bamStrand <- c(bamStrand,as.character(bam[[ii]]$strand))
        }
        return(bamStrand)
    }
getBamFull <-
    function(bamFile,rngs){
        para <- ScanBamParam(which=rngs,what=scanBamWhat(),reverseComplement=TRUE)
        bam <- scanBam(bamFile, param=para)
        return(bam)
    }
getBamSeqs <-
    function(bam){
        bamseqs <- DNAStringSet()
        for(ii in 1:length(bam)){
            bamseqs <- c(bamseqs,bam[[ii]]$seq)
        }
        return(bamseqs)
    }
getFasta <-
    function(fastaFile){
        fasta <- read.DNAStringSet(file=fastaFile,format='FASTA')
        return(fasta)
    }

trimFunc <-
    function(aligno,doMismatches){
        if(doMismatches){
            toDSS <- DNAStringSet(gsub('[ACGT]','M',aligno))
            mmt <- mismatchTable(aligno)
            ## for debugging
            if(sum(mmt$SubjectStart-mmt$SubjectEnd)>0){readline('mismatch table giving unexpected len')}
            ## end debuggin
            sp <- split(mmt,mmt$PatternId)
            namessp <- (as.numeric(names(sp)))

            for(i in seq_along(toDSS)){
                md <- Biostrings::as.matrix(toDSS[i])
                md[md %in% c('A','G','C','T')] <- 'M'
                if(i %in% namessp){
                    spsee <- sp[names(sp)==i]

                    md[spsee[[1]]$SubjectStart] <- as.character(spsee[[1]]$PatternSubstring)
                }
                toDSS[i] <- DNAStringSet(paste(md,collapse=''))

            }
            aligno <- toDSS
        }

        trim <- DNAStringSet(gsub('^-+','',aligno,perl=TRUE))
        trim <- DNAStringSet(gsub('-+$','',trim,perl=TRUE))
        return(trim)
    }
withflankNs <-
    function(aligno,trim){
        replb <- regexpr('^-+',aligno,perl=TRUE)
        reple <- regexpr('-+$',aligno,perl=TRUE)
        Nb <- DNAStringSet(DNAStringSet(rep(paste(rep("N",max(attr(replb,"match.length"))),collapse=''),length(replb))),start=1,width=sapply(attr(replb,"match.length"),function(x){max(0,x)} ))
        Ne <- DNAStringSet(DNAStringSet(rep(paste(rep("N",max(attr(reple,"match.length"))),collapse=''),length(reple))),start=1,width=sapply(attr(reple,"match.length"),function(x){max(0,x)}))
        return(DNAStringSet(paste(Nb,trim,Ne,sep='')))
    }

## add findSplitReadsOnly arg
mainAlignViewFull <-    function(bamFile,rngsRead,rngsAlign,filtSings=TRUE,findSplitReads=FALSE,filterbyMM=TRUE,MM=6,filterbyname=FALSE,filternames='',dedup=TRUE,typeArg ="global-local",substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],gapOpeningArg = -4, gapExtensionArg = -1,indelRate,mmRate,readLength,pairlimit=2e3,bamFileSplits=NULL,MMsplits=15,didSplits=FALSE,genomeName,findSplitReadsOnly=FALSE,doingContig=FALSE,verbose=FALSE){
    dftailA<-dftailB<-bamnamesF<-bamnamesF2<-bamnamesR<-bamnamesR2<-alignF<-alignF2<-alignR<-alignR2<-refInd<-bamnamesUse2<-bamnamesUseSplits<-refInd2<-alignedformatted<-list()
for(looptimes in 1:ifelse(doingContig,2,1)){
    if(doingContig) {rngsAlignUse<-rngsAlign[looptimes]} else rngsAlignUse<-rngsAlign
    ##print(doingContig);print(looptimes);print(rngsAlignUse)
### if findSplitReadsOnly and didSplits added
    if(findSplitReadsOnly==FALSE){
     alignedformatted <- doAlignAndformatFull(bamFile=bamFile,rngsRead=rngsRead,rngsAlign=rngsAlignUse,filtSings=filtSings,findSplitReads=findSplitReads,filterbyMM=filterbyMM,MM=MM,filterbyname=filterbyname,filternames=filternames,dedup=dedup,typeArg=typeArg,substitutionMat=substitutionMat,pairlimit=pairlimit,
                                gapOpeningArg=gapOpeningArg ,gapExtensionArg=gapExtensionArg,indelRate=indelRate,mmRate=mmRate,readLength=readLength,bamnamesRemove='',genomeName=genomeName,doingSplits=FALSE,pairedo='',findSplitReadsOnly=findSplitReadsOnly,didSplits=didSplits,verbose=verbose)
     if(length(alignedformatted)>0){


    dftailA[[looptimes]] <- alignedformatted[[1]]
    bamnamesF[[looptimes]]<-alignedformatted[[2]]
    bamnamesR[[looptimes]]<-alignedformatted[[3]]
    alignF[[looptimes]]<-alignedformatted[[4]]
    alignR[[looptimes]]<-alignedformatted[[5]]
     refInd[[looptimes]]<-alignedformatted[[6]]
        bamnamesUse2[[looptimes]] <- as.character(head(unique(dftailA[[looptimes]]$refname),-1))
        Fkeep2 <- which(bamnamesF[[looptimes]] %in% bamnamesUse2[[looptimes]])
        Rkeep2 <- which(bamnamesR[[looptimes]] %in% bamnamesUse2[[looptimes]])
        alignF[[looptimes]] <- alignF[[looptimes]][Fkeep2]
        alignR[[looptimes]] <- alignR[[looptimes]][Rkeep2]
        bamnamesF[[looptimes]] <- bamnamesF[[looptimes]][Fkeep2]
        bamnamesR[[looptimes]] <- bamnamesR[[looptimes]][Rkeep2]
    ## added for findSplitReadsOnly but now general purpose
    if(length(bamnamesF[[looptimes]])>0) names(bamnamesF[[looptimes]]) <- rep("F",length(bamnamesF[[looptimes]]))
    if(length(bamnamesR[[looptimes]])>0) names(bamnamesR[[looptimes]]) <- rep("R",length(bamnamesR[[looptimes]]))
    names(bamnamesUse2[[looptimes]])[bamnamesUse2[[looptimes]] %in% bamnamesF[[looptimes]]]='F'
    names(bamnamesUse2[[looptimes]])[bamnamesUse2[[looptimes]] %in% bamnamesR[[looptimes]]]='R'

    ## end added for findSplitReadsOnly


}else{

    dftailA[[looptimes]] <-data.frame()
    pval<-rep(1,readLength)
    bamnamesUse2[[looptimes]] <- vector()
    bamnamesF[[looptimes]] <-bamnamesR[[looptimes]] <-''
}
 }else{ ### else for findSplitReadsOnly
     bamnamesF[[looptimes]]<-bamnamesR[[looptimes]]<-''; dftailA[[looptimes]] <-data.frame();pval<-rep(1,readLength);bamnamesUse2[[looptimes]]<-vector();
 } ### added findSplitReadsOnly to if
    if(findSplitReadsOnly==TRUE | (length(unique(dftailA[[looptimes]]$ypostail))<5 & !is.null(bamFileSplits) & findSplitReads==TRUE) | didSplits==TRUE){
        ##print(rngsAlign);

        substitutionMatSplits<-nucleotideSubstitutionMatrix(match = 5, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)]

        pairedo<- ifelse(is.na(rownames(dftailA[[looptimes]])[1]),'',rownames(dftailA[[looptimes]])[1])
### added findSplitReadsOnly and didSplits

         alignedformattedSplits <-doAlignAndformatFull(bamFileSplits,rngsRead,rngsAlignUse,filtSings=filtSings,findSplitReads=findSplitReads,filterbyMM=filterbyMM,MM=MMsplits,filterbyname=filterbyname,filternames=filternames,dedup=dedup,typeArg=typeArg,substitutionMat=substitutionMatSplits,
                                pairlimit=pairlimit,gapOpeningArg=-100 ,gapExtensionArg= 0,indelRate=indelRate,mmRate=mmRate,readLength=readLength,
                                                       bamnamesRemove=c(bamnamesF[[looptimes]],bamnamesR[[looptimes]]),genomeName=genomeName,doingSplits=TRUE,pairedo=pairedo,findSplitReadsOnly=findSplitReadsOnly,didSplits=didSplits,verbose=verbose)
### added for findSplitReadsOnly
        if(length(alignedformattedSplits)==0 & findSplitReadsOnly){
         alignedformattedSplits <-doAlignAndformatFull(bamFile,rngsRead,rngsAlignUse,filtSings=filtSings,findSplitReads=findSplitReads,filterbyMM=filterbyMM,MM=MMsplits,filterbyname=filterbyname,filternames=filternames,dedup=dedup,typeArg=typeArg,substitutionMat=substitutionMatSplits,
                                pairlimit=pairlimit,gapOpeningArg=-100 ,gapExtensionArg= 0,indelRate=indelRate,mmRate=mmRate,readLength=readLength,
                                                       bamnamesRemove=c(bamnamesF[[looptimes]],bamnamesR[[looptimes]]),genomeName=genomeName,doingSplits=TRUE,pairedo=pairedo,findSplitReadsOnly=findSplitReadsOnly,didSplits=didSplits,verbose=verbose)
     }
### end added for findSplitReadsOnly



        if(length(alignedformattedSplits)>0){

        dftailB[[looptimes]]<-alignedformattedSplits[[1]]
        bamnamesF2[[looptimes]]<-alignedformattedSplits[[2]]
        bamnamesR2[[looptimes]]<-alignedformattedSplits[[3]]
        alignF2[[looptimes]]<-alignedformattedSplits[[4]]
        alignR2[[looptimes]]<-alignedformattedSplits[[5]]
        refInd2[[looptimes]]<-alignedformattedSplits[[6]]
        bamnamesUseSplits[[looptimes]]<- as.character(head(unique(dftailB[[looptimes]]$refname),-1))
        Fkeep2<-which(bamnamesF2[[looptimes]] %in% bamnamesUseSplits[[looptimes]])
        Rkeep2<-which(bamnamesR2[[looptimes]] %in% bamnamesUseSplits[[looptimes]])
        alignF2[[looptimes]]<- alignF2[[looptimes]][Fkeep2]
        alignR2[[looptimes]]<-alignR2[[looptimes]][Rkeep2]
        bamnamesF2[[looptimes]]<-bamnamesF2[[looptimes]][Fkeep2]
        bamnamesR2[[looptimes]]<-bamnamesR2[[looptimes]][Rkeep2]
        if(length(rngsAlignUse)==2){
        bestSplitScoreInd<-which.max(c(score(alignF2[[looptimes]]),score(alignR2[[looptimes]])))
        splitUse<-(c(DNAStringSet(as.character(pattern(alignF2[[looptimes]]))),DNAStringSet(as.character(pattern(alignR2[[looptimes]])))))[bestSplitScoreInd]
        splitUse<-gsub('-','',as.character(splitUse))
        names(splitUse)<-c(bamnamesF2[[looptimes]],bamnamesR2[[looptimes]])[bestSplitScoreInd]
        if(length(splitUse)==0) splitUse<-''
    }else splitUse<-''

        #pval2<-getpvalFull(bamnamesF2,bamnamesR2,alignF2,alignR2,readLength,refInd2,mmRate,indelRate)
        #pval<-pval*pval2
        bamnamesUse2[[looptimes]]<-c(bamnamesUse2[[looptimes]],bamnamesUseSplits[[looptimes]])
        ## added code for findSplitReadsOnly but generally applicable
        names(bamnamesUse2[[looptimes]])[bamnamesUse2[[looptimes]] %in% c(bamnamesF[[looptimes]],bamnamesF2[[looptimes]])]='F'
        names(bamnamesUse2[[looptimes]])[bamnamesUse2[[looptimes]] %in% c(bamnamesR[[looptimes]],bamnamesR2[[looptimes]])]='R'
        ## added code for findSplitReadsOnly but generally applicable
        if(length(bamnamesUse2[[looptimes]])==0) bamnamesUse2[[looptimes]] <- ''

        dftailB[[looptimes]]$ypostail[dftailB[[looptimes]]$ypostail != -1]=dftailB[[looptimes]]$ypostail[dftailB[[looptimes]]$ypostail != -1]+
            suppressWarnings(max(c(0,max(dftailA[[looptimes]]$ypostail))))
        dftailA[[looptimes]] <-unique(rbind(dftailA[[looptimes]],dftailB[[looptimes]]))
    }
#        dfs<-split(dftailA,dftailA$ypostail)


        didSplits<-TRUE
    }
}
    if(length(alignedformatted)>0) {
        if(doingContig) pval<-getpvalFull(bamnamesF[[1]],bamnamesR[[1]],alignF[[1]],alignR[[1]],readLength,refInd[[1]],mmRate,indelRate,bamnamesF[[2]],bamnamesR[[2]],alignF[[2]],alignR[[2]],refInd[[2]],findSplitReadsOnly)
        else pval<-getpvalFull(bamnamesF[[1]],bamnamesR[[1]],alignF[[1]],alignR[[1]],readLength,refInd[[1]],mmRate,indelRate,findSplitReadsOnly=findSplitReadsOnly)

    }
     if(didSplits){
         if(length(alignedformattedSplits)>0){
         if(doingContig) pval2<-getpvalFull(bamnamesF2[[1]],bamnamesR2[[1]],alignF2[[1]],alignR2[[1]],readLength,refInd2[[1]],mmRate,indelRate,bamnamesF2[[2]],bamnamesR2[[2]],alignF2[[2]],alignR2[[2]],refInd2[[2]],findSplitReadsOnly)
         else pval2<-getpvalFull(bamnamesF2[[1]],bamnamesR2[[1]],alignF2[[1]],alignR2[[1]],readLength,refInd2[[1]],mmRate,indelRate,findSplitReadsOnly=findSplitReadsOnly)

         pval<-pval*pval2
     }
     }

    if(exists('splitUse')==FALSE) splitUse<-''

result <- list(bamnames=bamnamesUse2[[1]],forplot=dftailA[[1]],pval=pval,didSplits=didSplits,splitRead=splitUse)
if(doingContig) result <-list(result, list(bamnames=bamnamesUse2[[2]],forplot=dftailA[[2]],pval=pval,didSplits=didSplits,splitRead=splitUse))

return(result)


}

getpvalFull<-function(bamnamesF,bamnamesR,alignF,alignR,readLength,refInd,mmRate,indelRate,bamnamesF2=character(),bamnamesR2=character(),
                      alignF2= PairwiseAlignments('A','A')[FALSE],alignR2=PairwiseAlignments('A','A')[FALSE],refInd2=matrix(),findSplitReadsOnly=FALSE){
        ## A little redundant

        forwardMtable <- table(mismatchTable(alignF)$PatternId)
                reverseMtable <- table(mismatchTable(alignR)$PatternId)

        forwardMtable2 <- table(mismatchTable(alignF2)$PatternId)
                reverseMtable2 <- table(mismatchTable(alignR2)$PatternId)

                allMismatchesF <- rep(0,length(alignF))
                allMismatchesR <- rep(0,length(alignR))

                allMismatchesF2 <- rep(0,length(alignF2))
                allMismatchesR2 <- rep(0,length(alignR2))


        allMismatchesF[as.numeric(names(forwardMtable))] <- forwardMtable
        names(allMismatchesF) <- rep('F',length(allMismatchesF))
        allMismatchesR[as.numeric(names(reverseMtable))] <- reverseMtable
        names(allMismatchesR)<-rep('R',length(allMismatchesR))


        allMismatchesF2[as.numeric(names(forwardMtable2))] <- forwardMtable2
        names(allMismatchesF2) <- rep('F2',length(allMismatchesF2))
        allMismatchesR2[as.numeric(names(reverseMtable2))] <- reverseMtable2
        names(allMismatchesR2)<-rep('R2',length(allMismatchesR2))

        if(findSplitReadsOnly){
                moreMism <- unlist(sapply(split(c(allMismatchesF,allMismatchesR,allMismatchesF2,allMismatchesR2),c(bamnamesF, bamnamesR,
                                                                                                                   bamnamesF2,bamnamesR2)),
                                          function(x){names(x)[which.min(x)]}))
            }else{
                moreMism <- unlist(sapply(split(c(allMismatchesF,allMismatchesR,allMismatchesF2,allMismatchesR2),c(bamnamesF, bamnamesR,
                                                                                                                   bamnamesF2,bamnamesR2)),
                                          function(x){x1s<-x[names(x) %in% c('F','R')];
                                                  x2s<-x[names(x) %in% c('F2','R2')];
                                                  x1max<-x1s[which.max(x1s)];
                                                  x2max<-x2s[which.max(x2s)];
                                                  xs<-c(x1max,x2max);
                                                  names(xs)[which.min(xs)]}))
            }
                
        alignFmm<-alignF[bamnamesF %in% names(moreMism)[moreMism=="F"]]
        bamnamesFmm<-bamnamesF[bamnamesF %in% names(moreMism)[moreMism=="F"]]

        alignF2mm<-alignF2[bamnamesF2 %in% names(moreMism)[moreMism=="F2"]]
        bamnamesF2mm<-bamnamesF2[bamnamesF2 %in% names(moreMism)[moreMism=="F2"]]


        tt<-table(start(subject(alignFmm)))
        tt2<-table(start(subject(alignF2mm)))
        if(any(tt>1)) {
            torem<-which(start(subject(alignFmm)) %in% as.numeric(names(tt))[tt>1])[cumsum(tt[tt>1])]
            alignFmm<-alignFmm[!(1:length(alignFmm) %in% torem)]
            bamnamesFmm<-bamnamesFmm[!(1:length(bamnamesFmm) %in% torem)]
        }
        if(any(tt2>1)) {
            torem2<-which(start(subject(alignF2mm)) %in% as.numeric(names(tt2))[tt2>1])[cumsum(tt2[tt2>1])]
            alignF2mm<-alignF2mm[!(1:length(alignF2mm) %in% torem2)]
            bamnamesF2mm<-bamnamesF2mm[!(1:length(bamnamesF2mm) %in% torem2)]
        }

        alignRmm<-alignR[bamnamesR %in% names(moreMism)[moreMism=="R"]]
        bamnamesRmm<-bamnamesR[bamnamesR %in% names(moreMism)[moreMism=="R"]]

        alignR2mm<-alignR2[bamnamesR2 %in% names(moreMism)[moreMism=="R2"]]
        bamnamesR2mm<-bamnamesR2[bamnamesR2 %in% names(moreMism)[moreMism=="R2"]]

        ttR<-table(start(subject(alignRmm)))
        ttR2<-table(start(subject(alignR2mm)))
        if(any(ttR>1)) {
            torem<-which(start(subject(alignRmm)) %in% as.numeric(names(ttR))[ttR>1])[cumsum(ttR[ttR>1])]
            alignRmm<-alignRmm[!(1:length(alignRmm) %in% torem)]
            bamnamesRmm<-bamnamesRmm[!(1:length(bamnamesRmm) %in% torem)]
        }
        if(any(ttR2>1)) {
            torem2<-which(start(subject(alignR2mm)) %in% as.numeric(names(ttR2))[ttR2>1])[cumsum(ttR2[ttR2>1])]
            alignR2mm<-alignR2mm[!(1:length(alignR2mm) %in% torem2)]
            bamnamesR2mm<-bamnamesR2mm[!(1:length(bamnamesR2mm) %in% torem2)]
        }


        ### new
        allMismatches1<-getMismatchFull(readLength,alignFmm,alignRmm)
        allMismatches2<-getMismatchFull(readLength,alignF2mm,alignR2mm)
        allIndel1<-getindelsFull(readLength,alignFmm,alignRmm,refInd)
        allIndel2<-getindelsFull(readLength,alignF2mm,alignR2mm,refInd2)
        ##allIndelRemInd<-allIndel[[1]]
        allIndel<-allIndel1[[1]]+allIndel2[[1]]
        allMismatches<-allMismatches1+allMismatches2
        ##browser()
        ###

        numreads<-length(alignFmm)+length(alignRmm)
        ##mmRate=rep(.01,100)
        mismatchBinom <- sapply(1:readLength,function(x){sum(dbinom(allMismatches[x]:numreads,numreads,mmRate[x]))})
        ##indelRate=rep(.005,100)

        indelBinom <- sapply(1:readLength,function(x){sum(dbinom(allIndel[x]:numreads,numreads,indelRate[x]))})

        pval <- indelBinom*mismatchBinom ##c(indelBinom,mismatchBinom)
        return(pval)


        result <- list(bamnamesUse2,dftailA,pval,allIndel,allMismatches,numreads,alignFmm,alignRmm,bamnamesFmm,bamnamesRmm)





}
## added findSplitReadsOnly argument
doAlignAndformatFull <-
    function(bamFile,rngsRead,rngsAlign,filtSings=TRUE,findSplitReads=FALSE,filterbyMM=TRUE,MM=6,filterbyname=FALSE,filternames='',dedup=TRUE,typeArg ="global-local",substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],
             pairlimit=2e3,gapOpeningArg = -4, gapExtensionArg = -1,indelRate,mmRate,readLength,bamnamesRemove='',genomeName='Hsapiens',doingSplits=FALSE,pairedo='',findSplitReadsOnly=FALSE,didSplits=FALSE,verbose=FALSE){



        bam <- getBamFull(bamFile,rngsRead)
        bamnames <- getBamNames(bam)
        bamseqs <- getBamSeqs(bam)

        if(length(bamseqs)==0) return(list())

        if(pairlimit != Inf){
            bampos <- getBamPos(bam)
            bamchr <- getBamChr(bam)
        }


        if(filternames[1] !=''){
            namekeepind<-which(bamnames %in% filternames)
            bamnames<-bamnames[namekeepind]
            bamseqs <- bamseqs[namekeepind]
            bampos <- bampos[namekeepind]
            bamchr <- bamchr[namekeepind]
        }

        if(dedup){
            dups <- which(Biostrings::duplicated(bamseqs))
            if(length(dups)>0){
                bamseqs <- bamseqs[-dups]
                bamnames <- bamnames[-dups]
                if(pairlimit != Inf){
                    bampos <- bampos[-dups]
                    bamchr <- bamchr[-dups]
                }
            }
        }


        if(bamnamesRemove[1] != ''){
            ##print(bamnamesRemove)

            indrem<-which(bamnames %in% bamnamesRemove)
            if(length(indrem>0)){
                bamseqs <-bamseqs[-indrem]
                bamnames <-bamnames[-indrem]
                bamchr<-bamchr[-indrem]
                bampos<-bampos[-indrem]
                ##print(indrem)
            }}

        if(length(bamseqs) > pairlimit){
            if(any(is.na(bampos))) bampos[which(is.na(bampos))]=301
            bamposR <- GRanges(seqnames=paste(bamchr,bamnames,sep=':'),ranges=IRanges(start=bampos-300,width=600))
            tt<-table(bamnames)
            nameskeep=names(tt)[tt==2]
            bamposR <- bamposR[!is.na(match(bamnames,nameskeep))]
            nameordR<-bamposR[order(as.character(gsub('(^.*?):','',seqnames(bamposR))))]

            tokeep<-which(width(pgap(ranges(nameordR[seq(1,length(nameordR),2)]),ranges(nameordR[seq(2,length(nameordR),2)])))>0)
            if(length(tokeep)==0) tokeep =1
            nameordRk<-c(nameordR[seq(1,length(nameordR),2)][tokeep],nameordR[seq(2,length(nameordR),2)][tokeep])
            tokeep2<-which(!is.na(match(bamnames,gsub('(^.*?:)', '',as.character(seqnames(nameordRk))))))
            bamseqs<-bamseqs[tokeep2]
            bamnames<-bamnames[tokeep2]
            bamseqs<-bamseqs[order(bamnames)][1:min(c(pairlimit,length(bamnames)))]
            bamnames<-bamnames[order(bamnames)][1:min(c(pairlimit,length(bamnames)))]
        }


        genome <- get(genomeName)
        if(length(rngsAlign)==1){
            refalign <- paste(suppressWarnings(Views(genome[[paste('chr',gsub('chr','',seqnames(rngsAlign)[1]),sep='')]],ranges(rngsAlign))),collapse='')
        }else if(length(rngsAlign)==2){
            rngsAlign<-sort(rngsAlign)

            ref1 <- suppressWarnings(Views(genome[[paste('chr',gsub('chr','',seqnames(rngsAlign)[1]),sep='')]],ranges(rngsAlign)[1]))
            ref2 <- suppressWarnings(Views(genome[[paste('chr',gsub('chr','',seqnames(rngsAlign)[2]),sep='')]],ranges(rngsAlign)[2]))
            refalign <- paste(ref1,ref2,sep='')
        }else{
            stop("More than 2 reference loci not supported")
        }



        refInd <- unlist(apply(cbind(start(rngsAlign),end(rngsAlign)),1,function(x){x[1]:x[2]}))
        refM <- unlist(apply(cbind(start(rngsAlign),end(rngsAlign)),1,function(x){rep(x[1],length(x[1]:x[2]))}))
        refNames <- unlist(apply(cbind(as.character(seqnames(rngsAlign)),width(rngsAlign)),1,function(x){rep(x[1],as.numeric(x[2]))}))
        refRa <- GRanges(seqnames=refNames,ranges=IRanges(start=refInd,width=1))

        alignFnew <- pairwiseAlignment(bamseqs,refalign,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)
        alignRnew <- pairwiseAlignment(reverseComplement(bamseqs),refalign,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)

        #### new for findSplitReadsOnly
        if(findSplitReadsOnly & didSplits==FALSE){
            delnamesF<-rep(1:length(alignFnew),lapply(deletion(alignFnew),length))
            delsF<-unlist(deletion(alignFnew))
            names(delsF)<-delnamesF
            startsF<-start(subject(alignFnew))[as.numeric(names(delsF))]
            splitsIndF<-as.numeric(names(delsF)[(refM[startsF+start(delsF)-1] != refM[startsF+end(delsF)-1])])
            delnamesR<-rep(1:length(alignRnew),lapply(deletion(alignRnew),length))
            delsR<-unlist(deletion(alignRnew))
            names(delsR)<-delnamesR
            startsR<-start(subject(alignRnew))[as.numeric(names(delsR))]
            splitsIndR<-as.numeric(names(delsR)[(refM[startsR+start(delsR)-1] != refM[startsR+end(delsR)-1])])
            alignF <- alignFnew[splitsIndF]
            alignR <- alignRnew[splitsIndR]
            bamnamesF <- bamnames[splitsIndF]
            bamnamesR <- bamnames[splitsIndR]

            #### end new for findSplitReadsOnly (added else statement below)
        }else{


        AlignBest <- apply(cbind(score(alignFnew),score(alignRnew)),1,which.max)

        alignF <- alignFnew[which(AlignBest==1)]
        alignR <- alignRnew[which(AlignBest==2)]

        bamnamesF <- bamnames[which(AlignBest==1)]
        bamnamesR <- bamnames[which(AlignBest==2)]
    } ### end else statement

        trimF <- trimFunc(alignF,TRUE)
        trimR <- trimFunc(alignR,TRUE)

        if(filterbyMM & findSplitReadsOnly==FALSE){
            doloop=TRUE
            MMuse<-MM

            while(doloop==TRUE){
            if(MMuse != MM & verbose) print(paste('mismatch limit increased to',MMuse,'to capture reads on both references'))
            if(MMuse > readLength){

                nogood<-which(refM[start(alignFnew@subject)] != refM[start(alignRnew@subject)])[1]
                if(is.na(nogood)) nogood=1
                alignFn<- alignFnew[nogood]
                alignRn<- alignRnew[nogood]
                bamnamesFn<-bamnames[nogood]
                bamnamesRn<-bamnames[nogood]
                trimFn <- trimFunc(alignFn,TRUE)
                trimRn <- trimFunc(alignRn,TRUE)
                Fkeep<-1;Rkeep<-1;bamnamesInit<-c(bamnamesFn,bamnamesRn);trimF<-trimFn;
                bamnamesF<-bamnamesFn;bamnamesR<-bamnamesRn;trimR<-trimRn;
                alignF<-alignFn;alignR<-alignRn;break
            }
            Fkeep <- which(rowSums(letterFrequency(trimF,c('A','C','T','G')))<MMuse)## & rowSums(letterFrequency(trimF,c('M')))>(.6*readLength))# | score(alignF)>0)
            Rkeep <- which(rowSums(letterFrequency(trimR,c('A','C','T','G')))<MMuse)## & rowSums(letterFrequency(trimR,c('M')))>(.6*readLength))# | score(alignR)>0)
            bamnamesInit <- (c(bamnamesF[Fkeep],bamnamesR[Rkeep]))
            if(length(rngsAlign)>1 &  bamnamesRemove[1]==''){
                if(length(bamnamesInit)>0){
                doloop<-all(aggregate(refM[c(start(alignF[Fkeep]@subject),start(alignR[Rkeep]@subject))],by=list(bamnamesInit),
                                    function(x){length(unique(x))})[,2]<2)
            }else{
                doloop<-TRUE
            }
            }else{
                doloop<-FALSE
            }
            MMuse<-MMuse+1
        }

            if(max(table(bamnamesInit))==1){
                Fkeep <- which(bamnamesF %in% bamnamesInit)
                Rkeep <- which(bamnamesR %in% bamnamesInit)
            }
        }else if(filterbyname){

            ## added code for findSplitReadsOnly
            Fkeep <- which(bamnamesF %in% filternames[!findSplitReadsOnly | names(filternames)=='F'])
            Rkeep <- which(bamnamesR %in% filternames[!findSplitReadsOnly | names(filternames)=='R'])
            ## end added code for findSplitReadsOnly
        }else if(findSplitReadsOnly){
            print('finding split reads only');
            Fkeep <- which(rowSums(letterFrequency(trimF,c('A','C','T','G')))<MM)## & rowSums(letterFrequency(trimF,c('M')))>(.6*readLength))# | score(alignF)>0)
            Rkeep <- which(rowSums(letterFrequency(trimR,c('A','C','T','G')))<MM)## & rowSums(letterFrequency(trimR,c('M')))>(.6*readLength))# | score(alignR)>0)
        }## end else for findSplitReadsOnly



        trimF <- trimF[Fkeep]
        trimR <- trimR[Rkeep]
        alignF <- alignF[Fkeep]
        alignR <- alignR[Rkeep]
        bamnamesF <- bamnamesF[Fkeep]
        bamnamesR <- bamnamesR[Rkeep]

        ### added for findSplitReadsOnly
        if(findSplitReadsOnly & (any(duplicated(bamnamesF)))){

            alignFord<-order(score(alignF),decreasing=TRUE)
            FkeepSplit<-alignFord[which(!duplicated(bamnamesF[alignFord]))]
            alignF<-alignF[FkeepSplit]
            trimF<-trimF[FkeepSplit]
            bamnamesF<-bamnamesF[FkeepSplit]

        }
        if(findSplitReadsOnly & (any(duplicated(bamnamesR)))){
            alignRord<-order(score(alignR),decreasing=TRUE)
            RkeepSplit<-alignRord[which(!duplicated(bamnamesR[alignRord]))]
            alignR<-alignR[RkeepSplit]
            trimR<-trimR[RkeepSplit]
            bamnamesR<-bamnamesR[RkeepSplit]
        }
        ## done for findSplitReadsOnly

        bamnamesUse <- c(bamnamesF,bamnamesR)##,bamnamesFsplit,bamnamesRsplit)

        correctPosFor <- start(alignF@subject)
        correctPosRev <- start(alignR@subject)

        deLF <- deletion(nindel(alignF))[,'WidthSum']
        deLR <- deletion(nindel(alignR))[,'WidthSum']



        frvec <- c(rep('F',length(trimF)),rep('R',length(trimR)))



        if((length(trimF)+length(trimR))==0) return(vector())

        dftailA <- dfmake(c(trimF,trimR),frvec  ,"fordf1A","bmat",c(correctPosFor,correctPosRev) ,c(deLF,deLR),bamnamesUse,
                          refalign,refInd,refRa,refM,filtSings,bamnamesRemove,findSplitReads,rngsAlign,doingSplits,pairedo)


        return(list(dftailA,bamnamesF,bamnamesR,alignF,alignR,refInd))
    }
