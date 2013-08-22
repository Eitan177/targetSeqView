\name{ViewAndScoreFull}
\alias{ViewAndScoreFull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ViewAndScoreFull
}
\description{
Perform full smith-waterman alignment for all three read-pair
configurations of interest (two supporting contiguous fragments, and one
supporting the structural variant) and return likelihood scores and a
data.frame for plotting.
}
\usage{
ViewAndScoreFull(filename,filtsings=TRUE,findSplitReads=FALSE,dedup=TRUE,allowedMM=6,initialExpansion=0,refexpansion=400,
         indelRate,mmRate,readLength,pairlimit=2e3,gapOpeningArg = -4,gapExtensionArg = -1,
         substitutionMat=nucleotideSubstitutionMatrix(match = 1,mismatch = -3)[c(1:4,8:9,15),c(1:4,8:9,15)],build='hg19',
         bsbuildprefix="BSgenome.Hsapiens.UCSC.",
         MMsplits=15,rngsAlign=GRanges(),validChr=c(1:22,'X','Y','M'),verbose=FALSE)

}
\arguments{
  \item{filename}{
the filename of the text file containing candidate SVs. Alternately, A data.frame or matrix with one row per candidate SV
}
  \item{initialExpansion}{
number of bases to add to the flanks of each loci before performing
SV-supporting alignment
}
\item{refexpansion}{
number of bases to add to the flanks of each loci before performing
contiguous-sequence-supporint alignment
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
\item{bamFilePath}{
full directory path of bam file (without the actual name of bam
file).
}
\item{allowedMM}{
the number of mismatches allowed in each read per read-pair aligned to the two
SV-supporting loci in order to be considered in
the score
}
\item{gapOpeningArg}{
the gap open penalty used in the pairwiseAlignment function call
}
\item{readLength}{
the length of reads used in the experiment
}
  \item{filtsings}{
boolean indicating if reads without a partner should be included in alignment/scoring
}
  \item{findSplitReads}{
boolean indicating if split-reads should be sought. Note, if this is set
to TRUE events must have a column indicating which bam file to use to
look for split-reads (named 'SplitsSample')
}
  \item{dedup}{
boolean indicating if reads with exactly matching sequences should be
considered only once (if TRUE) or as many times as they appear in the
bam file (if FALSE)
}
 \item{eventsToEval}{
a vector of rows to evaluate corresponding to the events (all rows are evaluated by default)
}
\item{MMsplits}{
the number of mismatches to allow in the alignment when searching for split-reads
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
A list containing two items per candidate; (1) a data.frame to pass to the
plotting function and (2) a likelihood score
}
\author{
Eitan Halper-Stromberg
}
\examples{
path <- system.file("extdata", package="targetSeqView")
nodes=1
registerDoMC(nodes)
candidateSVs<-new('candidates')
bamFilePath(candidateSVs) <- path
candidatesFileName(candidateSVs) <-
file.path(path,'targetCaptureSVs.txt')
build(candidateSVs) <- 'hg19'
readLength(candidateSVs) <- 100
mmRate(candidateSVs) <- precomputedTargetCapture100bpMMRate()
indelRate(candidateSVs) <- precomputedTargetCapture100bpIndelRate()

     retvals <-ViewAndScoreFull(filename=candidateSVs@candidatesFileName,
                                  readLength=candidateSVs@readLength,
                                  bamFilePath=candidateSVs@bamFilePath,
                                  mmRate=candidateSVs@mmRate,
                                  indelRate=candidateSVs@indelRate,
eventsToEval=2,
,verbose=TRUE)
     candidateSVs@forplot <- retvals[['forplot']]
     candidateSVs@fullScore <- retvals[['score']]
 ## this function is called within the fullScoreAndView method
}
