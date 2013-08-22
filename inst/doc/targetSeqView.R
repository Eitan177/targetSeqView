### R code from vignette source 'targetSeqView.Rnw'

###################################################
### code chunk number 1: <setup
###################################################
    library(targetSeqView)
    options(width=70)


###################################################
### code chunk number 2: pkgs
###################################################



###################################################
### code chunk number 3: instantiateDel
###################################################
path <- system.file("extdata", package="targetSeqView")

## This method utilizes the foreach package for parallelization, set nodes to however many cpus are
## available.
nodes=1
registerDoMC(nodes)

## create an instance of the candidates class
candidateDels<-new('candidates')
## set the path where bam files are located (if not in the currect working directory)
bamFilePath(candidateDels)<-path
## set the name of the text file containing candidate SVs (full path if not in the working directory)
candidatesFileName(candidateDels)<-file.path(path,'wholeGenomeDeletionCandidates.txt')
## set the build of the (human) genome
build(candidateDels) <- 'hg19'
## set the read length
readLength(candidateDels) <- 101
## set the mismatch rate for each position along the read length
mmRate(candidateDels) <- precomputedWholeGenome101bpMMRate()
## set the indel rate for reach position along the read length
indelRate(candidateDels) <- precomputedWholeGenome101bpIndelRate()


###################################################
### code chunk number 4: computeErrorRates (eval = FALSE)
###################################################
##  normalBam <- 'Path/To/Normal/bamfile.bam'
## errorRates<-getErrorRate(normalBam)
## mmRate(candidateDels) <- errorRates[['mmRate']]
## indelRate(candidateDels) <- errorRates[['indelRate']]


###################################################
### code chunk number 5: quickscore
###################################################
candidateDels<- quickScore(candidateDels,verbose=TRUE)
## view values returned
print(candidateDels@quickScore)
### In this case we have validation data for these candidates
indexOfvalidated <-1:10
validated<-candidateDels@quickScore[indexOfvalidated]

failedvalidation<-candidateDels@quickScore[-indexOfvalidated]


###################################################
### code chunk number 6: plotquickscore
###################################################
boxplot(list(validated=validated,failed=failedvalidation,
    all=candidateDels@quickScore),ylab='log likelihood score')


###################################################
### code chunk number 7: instantiateSV
###################################################
## create an instance of the candidates class
candidateSVs<-new('candidates')
bamFilePath(candidateSVs) <- path
candidatesFileName(candidateSVs) <- file.path(path,'targetCaptureSVs.txt')
build(candidateSVs) <- 'hg19'
readLength(candidateSVs) <- 100
mmRate(candidateSVs) <- precomputedTargetCapture100bpMMRate()
indelRate(candidateSVs) <- precomputedTargetCapture100bpIndelRate()

## fullScoreAndview will perform full smith-waterman realignment for all reads in the 3
## configurations. In addition, if the input text file contains a SplitsSample column,
## the function will look for split-reads within the bam file specified by the column 'SplitsSample'
candidateSVs<-fullScoreAndView(candidateSVs,verbose=TRUE,findSplitReads=TRUE)


###################################################
### code chunk number 8: printscore
###################################################
print(candidateSVs@fullScore)


###################################################
### code chunk number 9: plot1
###################################################
plotSV(candidateSVs,indices=1,flipLeftandRight=TRUE,pdfname='fig1.pdf')


###################################################
### code chunk number 10: plot2
###################################################
plotSV(candidateSVs,indices=2,pdfname='fig2.pdf')


###################################################
### code chunk number 11: plot3
###################################################
plotSV(candidateSVs,indices=3,flipLeftandRight=TRUE,pdfname='fig3.pdf',width=10)


###################################################
### code chunk number 12: sessionInfo
###################################################
toLatex(sessionInfo())


