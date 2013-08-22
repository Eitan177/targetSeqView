\name{doQuickScore}
\alias{doQuickScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
doQuickScore
}
\description{
Perform smith-waterman alignment for the two read-pair
configurations supporting contiguous fragments and use the cigar
strings to infer alignments for the configuration supporting the SV.
Return likelihood scores only.
}
\usage{
doQuickScore(filename,
             readLength,bamFilePath='',gapOpeningArg = -4,
             dedup=TRUE,eventsToEval=Inf,mmRate,indelRate,pairlimit=2e3,
             gapExtensionArg = -1,substitutionMat=nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],
             rngsAlign=GRanges(),build='hg19',bsbuildprefix="BSgenome.Hsapiens.UCSC.",conservativeContigAlign=FALSE,validChr=c(1:22,'X','Y','M'),verbose=FALSE)
}
\arguments{
  \item{filename}{
the filename of the text file containing candidate SVs. Alternately, A data.frame or matrix with one row per candidate SV
}
\item{readLength}{
the length of reads used in the experiment
}
\item{bamFilePath}{
full directory path of bam file (without the actual name of bam
file).
}
\item{gapOpeningArg}{
the gap open penalty used in the pairwiseAlignment function call
}
  \item{dedup}{
boolean indicating if reads with exactly matching sequences should be
considered only once (if TRUE) or as many times as they appear in the
bam file (if FALSE)
}
 \item{eventsToEval}{
a vector of rows to evaluate corresponding to the events (all rows are evaluated by default)
}
\item{mmRate}{
a vector holding the position-specific mismatch rate in your experiment
(which may be calculated by calling 'getErrorRate'). The length of this
vector should be the same as the length of reads used in the experiment
}
\item{indelRate}{
a vector holding the position-specific indel rate in your experiment
(which may be calculated by calling 'getErrorRate'). The length of this
vector should be the same as the length of reads used in the experiment
}
 \item{pairlimit}{
the maximum number of read-pairs to consider for alignment/scoring.
}
\item{gapExtensionArg}{
the gap extension penalty used in the pairwiseAlignment function call
}
\item{substitutionMat}{
the nucleotideSubstitutionMatrix passed to the pairwiseAlignment
function call.
}
\item{rngsAlign}{
list of GRanges objects to use for alignment in place of loci read in
from text file
}
\item{build}{
genome build (this should be the same as what was used to align reads in
the bam file)
}
\item{bsbuildprefix}{
the prefix of the biostrings reference object
}
\item{validChr}{
a vector of the names of the chromosomes from the species of interest
(unimplemented here)
}
\item{verbose}{
print out progress messages
}
}
\details{
see example
}
\value{
A vector of likelihood scores, one for each evaluated candidate SV
}
\author{
Eitan Halper-Stromberg
}
\examples{
path <- system.file("extdata", package="targetSeqView")
nodes=1
registerDoMC(nodes)

## create an instance of the candidates class
candidateDels<-new('candidates')
## set the path where bam files are located (if not in the currect
working directory)
bamFilePath(candidateDels)<-path
## set the name of the text file containing candidate SVs (full path if
not in the working directory)
candidatesFileName(candidateDels)<-file.path(path,'wholeGenomeDeletionCandidates.txt')
## set the build of the (human) genome
build(candidateDels) <- 'hg19'
## set the read length
readLength(candidateDels) <- 101
## set the mismatch rate for each position along the read length
mmRate(candidateDels) <- precomputedWholeGenome101bpMMRate()
## set the indel rate for reach position along the read length
indelRate(candidateDels) <- precomputedWholeGenome101bpIndelRate()

    candidateDels@quickScore<-doQuickScore(filename=candidateDels@candidatesFileName,
                                  readLength=candidateDels@readLength,
                                  bamFilePath=candidateDels@bamFilePath,
                                  mmRate=candidateDels@mmRate,
                                  indelRate=candidateDels@indelRate,
                                  eventsToEval=1:2,verbose=TRUE)
 ## this function is called within the quickScore method
}
