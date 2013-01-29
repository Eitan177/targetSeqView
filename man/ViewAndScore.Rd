\name{ViewAndScore}
\alias{ViewAndScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ViewAndScore
}
\description{
get reads from a bam file, realign reads in three configurations, and return an object with realignments and
likelihood scores
}
\usage{
ViewAndScore(filename, initialExpansion = 500, estimateIndelRate = TRUE, indelRate = 0.005, estimateMmRate = TRUE, mmRate = 0.01, getReadLength = TRUE, readLength = 100, scoreOnly = FALSE, normalBam = "", bamFilePath = "", build="hg19",verbose = TRUE)
}
\arguments{
  \item{filename}{
    The bam file holding alignments to read in
  }
  \item{initialExpansion}{
    expand the reference sequence by this number of nucleotides on both the
    5' and 3' ends, relative to the coordinates read in from the text file
    before realigning reads in the configuraiton supporting the structural variant
  }
  \item{estimateIndelRate}{
    if true estimate the indel rate from good alignments taken from a bam
    file with normal alignments
  }
  \item{indelRate}{
    input indel rate if this is not to be estimated from normal alignments
    in a supplied bam file
  }
  \item{estimateMmRate}{
    if true estimate the mismatch rate from good alignments taken from a bam
    file with normal alignments
  }
  \item{mmRate}{
    input mismatch rate if this is not to be estimated from normal alignments
    in a supplied bam file
  }
  \item{getReadLength}{
    if true get the readLength from the bam file with normal alignments
  }
  \item{readLength}{
    input read length instead of getting it from a bam file with normal alignments
  }
  \item{scoreOnly}{
    depreciated
  }
  \item{normalBam}{
    name of the normal bam file (with full path) to use for estimating
    the mismatch and indel rates
  }
  \item{bamFilePath}{
    path to bam file with reads to realign
  }
  \item{build}{
    human genome referene build to use
  }
  \item{verbose}{
    output progress
  }
}
\details{
see example
}
\value{
  \item{list element 1 }{object holding information necessary to plot realignments}
  \item{list element 2 }{likelihood scores}

}
\author{
Eitan Halper-Stromberg
}
\examples{
library(grid)
path <- system.file("extdata", package="targetSeqView")
filename=file.path(path, "twoSVJunctionsFailed.txt")
nodes=1
registerDoMC(nodes)
retfail=ViewAndScore(filename=filename,bamFilePath=path,
estimateIndelRate=FALSE,estimateMmRate=FALSE,getReadLength=FALSE,
build='hg19',verbose=TRUE)
print('The likelihood scores for the events failing validation:')
print(retfail[[2]])
p1=formatPlot(retfail[[1]][[1]][[1]][[2]],title='Alignment supporting a structural variant')
p2=formatPlot(retfail[[1]][[1]][[2]][[2]],title='Alignment supporting no structural variant')
p3=formatPlot(retfail[[1]][[1]][[3]][[2]],title='Alignment supporting no structural variant')
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 1)))
print(p1,vp = viewport(layout.pos.row = 1, layout.pos.col=1))
print(p2,vp = viewport(layout.pos.row = 2, layout.pos.col=1))
print(p3,vp = viewport(layout.pos.row = 3, layout.pos.col=1))

filename=file.path(path,"twoSVJunctionsPassed.txt")
retpass=ViewAndScore(filename=filename,bamFilePath=path,
estimateIndelRate=FALSE,estimateMmRate=FALSE,getReadLength=FALSE,
build='hg19',verbose=TRUE)
print('The likelihood scores for the events passing validation:')
print(retpass[[2]])
p1=formatPlot(retpass[[1]][[1]][[1]][[2]],title='Alignment supporting a structural variant')
p2=formatPlot(retpass[[1]][[1]][[2]][[2]],title='Alignment supporting no structural variant')
p3=formatPlot(retpass[[1]][[1]][[3]][[2]],title='Alignment supporting no structural variant')
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 1)))
print(p1,vp = viewport(layout.pos.row = 1, layout.pos.col=1))
print(p2,vp = viewport(layout.pos.row = 2, layout.pos.col=1))
print(p3,vp = viewport(layout.pos.row = 3, layout.pos.col=1))

}
