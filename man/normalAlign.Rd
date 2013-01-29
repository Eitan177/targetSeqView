\name{normalAlign}
\alias{normalAlign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
normalAlign
}
\description{
realign normal readpairs from a bam file. Normal readpairs are defined
as readspairs that align in the proper orientation and the proper
fragment length as defined by the flags in the bam file
}
\usage{
normalAlign(bamFile, typeArg = "global-local", substitutionMat =
nucleotideSubstitutionMatrix(match = 1, mismatch = -3)[c(1:4, 8:9, 15),
c(1:4, 8:9, 15)], gapOpeningArg = -4, gapExtensionArg = -1,build = "hg19")
}
\arguments{
  \item{bamFile}{
name of bam file with normal alignments
}
  \item{typeArg}{
method for realignment, passed to pairwiseAlignment function
}
  \item{substitutionMat}{
substitution matrix, passed to pairwiseAlignment function
}
  \item{gapOpeningArg}{
penalty for opening a gap, passed to pairwiseAlignment function
}
  \item{gapExtensionArg}{
penalty to extend a gap, passed to pairwiseAlignment function
}
  \item{build}{
human genome referene build to use
}
}

\details{
see example
}
\value{
PairwiseAlignmentsSingleSubject  
}
\author{
Eitan Halper-Stromberg
}
\examples{
path <- system.file("extdata", package="targetSeqView")
normalBam=file.path(path, "tmp.bam")
    norms=normalAlign(bamFile=normalBam)
    rL=max(width((DNAStringSet(gsub('-','',as.character(pattern(norms)))))))
    readLength=rL
    print(readLength)
    indls=nindel(norms)
    indelRateEstimate=mean(indls@insertion[,2]+indls@deletion[,1])/readLength
    indelRate=indelRateEstimate
    print(indelRate)
    mms=rep(0,length(norms))
    mms[as.numeric(names(table(mismatchTable(norms)[,1])))]=table(mismatchTable(norms)[,1])
    mmRateEstimate=mean(mms)/readLength
    mmRate=mmRateEstimate
    print(mmRate)
}