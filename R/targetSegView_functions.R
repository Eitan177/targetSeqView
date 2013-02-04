
dfmake2 <-
    function(fordf1A,deL,bamn,refalign,refInd,refM,filtSings,rngsAlign){
        ## combine paired end

        fordPEs <- split(fordf1A,bamn)

        if(filtSings){
            ## make sure two reads pass per read name pass to this point
            btab <- table(bamn);fordPEs <- fordPEs[names(fordPEs) %in% names(btab)[btab>1]]
            ## make sure for each pair of read, majority of reference is on different reference if there are two reference
            if(length(rngsAlign)>1){
                difrefs <- sapply(fordPEs,function(x){freadref=table(x[[1]][,3]);sreadref=table(x[[2]][,3]);names(freadref)[which.max(freadref)] != names(sreadref)[which.max(sreadref)]})
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

        fordPEsCombDup <- lapply(fordPEsComb,function(x){x2=x[order(as.numeric(x[,2]),x[,1]),];x3=x2[!duplicated(x2[,2]),];if(length(which(duplicated(x2[,2])))>(nrow(x3)*.7) & length(rngsAlign)>1){NA}else{x3}})
        fordPEsCombDup[which(is.na(fordPEsCombDup))]<-NULL
        vvA <- unlist(lapply(fordPEsCombDup,function(x){as.numeric(as.character(data.frame(x)[1,2]))}))


        fordf2A <- sapply(order(as.numeric(vvA)),function(x){list(fordPEsCombDup[[x]])})
        fordf3A <- sapply(1:length(fordf2A),function(x){list(cbind(fordf2A[[x]],rep((x),nrow(fordf2A[[x]]))))})
        names(fordf3A) <- names(fordPEsCombDup)[(order(as.numeric(vvA)))]
        namecol <- rep(names(fordf3A),unlist(sapply(fordf3A,nrow)))
        dftailA <- data.frame(do.call(rbind,fordf3A),namecol)
        for(ii in 2:4){
            dftailA[,ii] <- as.numeric(as.character(dftailA[,ii]))
        }

        colnames(dftailA) <- c('afvec2','posvec2','ref','ypostail','refname')

        prerefrange <- aggregate(dftailA$posvec2,by=list(dftailA$ref),FUN=range)
        prerefsee <- substring(refalign,match(prerefrange$x[,1],refInd),match(prerefrange$x[,2],refInd))
        indr <- eval(parse(text=paste('c(',paste(apply(cbind(match(prerefrange$x[,1],refInd),match(prerefrange$x[,2],refInd)),1,function(x){paste(x,collapse=':')}),collapse=','),')')))
        refIndsee <- refInd[indr]
        refSide <- refM[indr]
        refsee <- Biostrings::unlist(DNAStringSet(prerefsee))

        formatrefsee <- data.frame(t(Biostrings::as.matrix(DNAStringSet(refsee))),refIndsee,refSide, -1,'refseq')
        colnames(formatrefsee) <- colnames(dftailA)
        dftailA <- rbind(dftailA,formatrefsee)#,zeros)

        return(dftailA)
    }
dfmake <-
    function(aa,frvec,fname,amatname,startpos,deL,bamn,refalign,refInd,refM,filtSings,rngsAlign){

        ees <- rep('E',7)
        fordf1A <- mapply(function(x,y,z){x2=strsplit(gsub('S','N',toString(x)),'')[[1]];if(z=='F'){newy=y; x2=c(x2,ees) }else{newy=max(1,y-7); x2=c(rep('E',y-newy),x2)};
                                          ## na.omit is to treat condition where ees run off the end
                                          cc=na.omit(cbind(x2,refInd[newy:(length(x2)+newy-1)],refM[newy:(length(x2)+newy-1)]));dd=diff(as.numeric(cc[,2]));
                                          ll=list(cc)},Biostrings::as.list(aa),as.list(startpos),as.list(frvec))
        return(dfmake2(fordf1A,deL,bamn,refalign,refInd,refM,filtSings,rngsAlign))
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
getBam <-
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
getRandomBamRecords <-
    function(bamFile,reverseC=FALSE){
        nGoodAlignments <- system(paste("echo `samtools flagstat", bamFile,"| grep 'properly paired' | awk '{print $1}'`"),intern=TRUE)
        fracAlignmentsGet <- system(paste('echo "scale=6; 10000/',nGoodAlignments,'" | bc',sep=''),intern=TRUE)
        system(paste('samtools view -hf 0x2',bamFile,'| samtools view -bhS -s',fracAlignmentsGet,'- > tmp.bam'))
        para <- ScanBamParam(what=scanBamWhat(),reverseComplement=reverseC)
        bam <- scanBam('tmp.bam',param=para)
        return(bam)
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

mainAlignView <-
    function(bamFile,rngsRead,rngsAlign,filtSings=TRUE,filterbyMM=TRUE,MM=6,filterbyname=FALSE,returnScoreOnly=FALSE,filternames,dedup=TRUE,typeArg ="global-local",substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],gapOpeningArg = -4, gapExtensionArg = -1,indelRate=.005,mmRate=.01,readLength=100){

        bam <- getBam(bamFile,rngsRead)
        bamnames <- getBamNames(bam)
        bamseqs <- getBamSeqs(bam)
        if(dedup){
            dups <- which(Biostrings::duplicated(bamseqs))
            if(length(dups)>0){
                bamseqs <- bamseqs[-dups]
                bamnames <- bamnames[-dups]
            }
        }


        rngsAlign <- rngsAlign#+100
        Hsapiens <- get("Hsapiens")
        if(length(GenomicRanges::unique(seqnames(rngsAlign)))==1){
            refalign <- paste(suppressWarnings(Views(Hsapiens[[paste('chr',seqnames(rngsAlign)[1],sep='')]],ranges(rngsAlign))),collapse='')
        }else if(length(GenomicRanges::unique(seqnames(rngsAlign)))==2){

            ref1 <- suppressWarnings(Views(Hsapiens[[paste('chr',seqnames(rngsAlign)[1],sep='')]],ranges(rngsAlign)[1]))
            ref2 <- suppressWarnings(Views(Hsapiens[[paste('chr',seqnames(rngsAlign)[2],sep='')]],ranges(rngsAlign)[2]))
            refalign <- paste(ref1,ref2,sep='')
        }else{
            stop("More than 2 reference loci not supported")
        }



        refInd <- unlist(apply(cbind(start(rngsAlign),end(rngsAlign)),1,function(x){x[1]:x[2]}))
        refM <- unlist(apply(cbind(start(rngsAlign),end(rngsAlign)),1,function(x){rep(x[1],length(x[1]:x[2]))}))



        mat <- substitutionMat;
        alignFnew <- pairwiseAlignment(bamseqs,refalign,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)
        alignRnew <- pairwiseAlignment(reverseComplement(bamseqs),refalign,type =typeArg,substitutionMatrix=substitutionMat,gapOpening = gapOpeningArg, gapExtension = gapExtensionArg)


        AlignBest <- apply(cbind(score(alignFnew),score(alignRnew)),1,which.max)

        alignF <- alignFnew[which(AlignBest==1)]
        alignR <- alignRnew[which(AlignBest==2)]

        bamnamesF <- bamnames[which(AlignBest==1)]
        bamnamesR <- bamnames[which(AlignBest==2)]



        trimF <- trimFunc(alignF,TRUE)
        trimR <- trimFunc(alignR,TRUE)



        if(filterbyMM){

            Fkeep <- which(rowSums(letterFrequency(trimF,c('A','C','T','G')))<MM)# | score(alignF)>0)
            Rkeep <- which(rowSums(letterFrequency(trimR,c('A','C','T','G')))<MM)# | score(alignR)>0)
            bamnamesInit <- (c(bamnamesF[Fkeep],bamnamesR[Rkeep]))
            if(max(table(bamnamesInit))==1){
                Fkeep <- which(bamnamesF %in% bamnamesInit)
                Rkeep <- which(bamnamesR %in% bamnamesInit)
            }
        }else if(filterbyname){
            Fkeep <- which(bamnamesF %in% filternames)
            Rkeep <- which(bamnamesR %in% filternames)
        }

        trimF <- trimF[Fkeep]
        trimR <- trimR[Rkeep]
        alignF <- alignF[Fkeep]
        alignR <- alignR[Rkeep]
        bamnamesF <- bamnamesF[Fkeep]
        bamnamesR <- bamnamesR[Rkeep]
        bamnamesUse <- c(bamnamesF,bamnamesR)

        correctPosFor <- start(alignF@subject)
        correctPosRev <- start(alignR@subject)
        deLF <- deletion(nindel(alignF))[,'WidthSum']
        deLR <- deletion(nindel(alignR))[,'WidthSum']



        frvec <- c(rep('F',length(trimF)),rep('R',length(trimR)))




        dftailA <- dfmake(c(trimF,trimR),frvec  ,"fordf1A","bmat",c(correctPosFor,correctPosRev) ,c(deLF,deLR),bamnamesUse,refalign,refInd,refM,filtSings,rngsAlign)


        ## exlude refseq, hence the head(..,-1)
        bamnamesUse2 <- as.character(head(unique(dftailA$refname),-1))
        Fkeep2 <- which(bamnamesF %in% bamnamesUse2)
        Rkeep2 <- which(bamnamesR %in% bamnamesUse2)

        ## A little redundant
        trimF <- trimF[Fkeep2]
        trimR <- trimR[Rkeep2]
        alignF <- alignF[Fkeep2]
        alignR <- alignR[Rkeep2]
        bamnamesF <- bamnamesF[Fkeep2]
        bamnamesR <- bamnamesR[Rkeep2]



        forwardMtable <- table(mismatchTable(alignF)$PatternId)
        reverseMtable <- table(mismatchTable(alignR)$PatternId)
        allMismatchesF <- rep(0,length(alignF))
        allMismatchesR <- rep(0,length(alignR))
        allMismatchesF[as.numeric(names(forwardMtable))] <- forwardMtable
        allMismatchesR[as.numeric(names(reverseMtable))] <- reverseMtable
        allMisM <- unlist(sapply(split(c(allMismatchesF,allMismatchesR),c(bamnamesF, bamnamesR)),max))
        indelF <- nindel(alignF)
        ## For insertions we want every base involved in the insertion because those are read bases
        ## For deletions we want only the lengths because bases are from the reference, length will tell us that a read base is involved
        allindelF <- indelF@insertion[,2]+indelF@deletion[,1]
        indelR <- nindel(alignR)
        allindelR <- indelR@insertion[,2]+indelR@deletion[,1]
        allindel <- unlist(sapply(split(c(allindelF,allindelR),c(bamnamesF, bamnamesR)),max))
        indelBinom <- sapply(allindel,function(x){max(sum(dbinom(x:readLength,readLength,indelRate)),  1-ifelse(x==0,0,sum(dbinom(0:(x-1),readLength,indelRate))))})
        mismatchBinom <- sapply(allMisM,function(x){max(sum(dbinom(x:readLength,readLength,indelRate)), 1-ifelse(x==0,0,sum(dbinom(0:(x-1),readLength,mmRate))))})
        pval <- indelBinom*mismatchBinom
        result <- list(bamnamesUse2,dftailA,pval,allindel,allMisM)

        return(result)
    }

