\name{alignView}
\alias{alignView}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
alignView
}
\description{
get reads from a bam file, realign reads in three configurations, and return an object with realignments and
alignment probabilities which will be used by the ViewAndScore function
to make likelihoods
}
\usage{
alignView(events, filtsings = TRUE, scoreOnly = FALSE, allowedMM = 6, initialExpansion = 0, refexpansion = 400, indelRate = 0.005, mmRate = 0.01, readLength = 100, build="hg19",verbose = FALSE)
}
\arguments{
  \item{events}{
    The bam file holding alignments to read in
  }
  \item{filtsings}{
    filter reads that map to one of the two loci but whos partner does not
  }
  \item{scoreOnly}{
    depreciated
  }
  \item{allowedMM}{
    filter reads aligning with more than this many mismatches to the
    rearranged loci 
  }
  \item{initialExpansion}{
    expand the reference sequence by this number of nucleotides on both the
    5' and 3' ends, relative to the coordinates read in from the text file
    before realigning reads in the configuraiton supporting the structural variant
  }
  \item{refexpansion}{
    expand the reference sequence by this number of nucleotides on both the
    5' and 3' ends, relative to the coordinates read in from the text file
    before realigning reads in the configurations not supporting the
    strucutural variant
  }
  \item{indelRate}{
    input indel rate
  }
  \item{mmRate}{
    input mismatch rate
  }
  \item{readLength}{
    read length
  }
  \item{build}{
    build of human genome reference to use
  }
  \item{verbose}{
    output progress
  }
}
\details{
see example
}
\value{
object holding one list element per candidate rearrangement
}
\author{
Eitan Halper-Stromberg
}
\examples{
path <- system.file("extdata", package="targetSeqView")
filename=file.path(path, "twoSVJunctionsFailed.txt")
events=as.matrix(read.delim(file=filename,header=TRUE))
events[,'Sample']=file.path(path,events[,'Sample'])
events[,'Chr1']=gsub('chrchr','chr',tolower(paste('chr',gsub('','',events[,'Chr1']),sep='')))
events[,'Chr2']=gsub('chrchr','chr',tolower(paste('chr',gsub('','',events[,'Chr2']),sep='')))
nodes=1
registerDoMC(nodes)
retfail=alignView(events=events,initialExpansion=500,build='hg19')
for(ii in 1:length(retfail)){
    dev.new();
    formatPlot(retfail[[ii]][[1]][[2]],title='Alignment supporting a structural variant')
    dev.new();
    formatPlot(retfail[[ii]][[2]][[2]],title='Alignment supporting no structural variant')
    dev.new();
    formatPlot(retfail[[ii]][[3]][[2]],title='Alignment supporting no structural variant')
}
}

