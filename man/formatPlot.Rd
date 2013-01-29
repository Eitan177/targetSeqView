\name{formatPlot}
\alias{formatPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
formatPlot
}
\description{
plot realigned reads
}
\usage{
formatPlot(dftailA,title)
}
\arguments{
  \item{dftailA}{
    data frame containing output from alignView or ViewAndScore
  }
  \item{title}{
    title of plot
  }
}
\details{
  see example
}
\value{
  a ggplot2 object
}
\author{
  Eitan Halper-Stromberg
}
\examples{
path <- system.file("extdata", package="targetSeqView")
filename=file.path(path, "twoSVJunctionsFailed.txt")
nodes=1
registerDoMC(nodes)
retfail=ViewAndScore(filename=filename,bamFilePath=path,
estimateIndelRate=FALSE,estimateMmRate=FALSE,getReadLength=FALSE,
build='hg19',verbose=TRUE)
print('The likelihood scores for the events failing validation:')
print(retfail[[2]])
for(ii in 1:length(retfail[[1]])){
    dev.new();
    formatPlot(retfail[[1]][[ii]][[1]][[2]],title='Alignment supporting a structural variant')
    dev.new();
    formatPlot(retfail[[1]][[ii]][[2]][[2]],title='Alignment supporting no structural variant')
    dev.new();
    formatPlot(retfail[[1]][[ii]][[3]][[2]],title='Alignment supporting no structural variant')
}
}