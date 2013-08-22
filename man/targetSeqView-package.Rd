\name{targetSeqView-package}
\alias{targetSeqView-package}
\alias{targetSeqView}
\docType{package}
\title{
targetSeqView-Realignment and likelihood scoring for candidate
structural variants
}
\description{
Score candidate structural variants based upon alignment likelihood. The
user must provide a text file with one row per candidate. Each row
should contain the start and end of the two loci involved in the event
and the name of the bam file containing aligned reads. The package loads
reads-aligning within the indicated loci and performs smith-waterman
realignment in 3 configurations, 2 supporting contiguous sequences and 1
supporting the structural variant. A probability is calcuated for each
alignment and the likelihood of an event is constructed based upon
weighing the probability that reads were generated from the structural
variants vs. the probability that reads were generated from a contiguous
sequence. Alignment views are also generated for manual inspection.
}
\details{
\tabular{ll}{
Package: \tab targerSeqView\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-08-21\cr
License: \tab Johsn Hopkins School of Medicine\cr
}
see the vignette for details
}
\author{
Eitan Halper-Stromberg
Maintainer: Eitan Halper-Stromberg
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
type ?ViewAndScoreFull.Rd at the command prompt for an example
}
