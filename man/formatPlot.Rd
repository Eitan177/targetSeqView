\name{formatPlot}
\alias{formatPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
formatPlot
}
\description{
create pictures for the three alignments considered in the construction
of a likelihood score.
}
\usage{
formatPlot(dfplot,mismatchcolor='red',title='',covMaxDisplay=Inf,flipLeftandRight=FALSE,cropGaps=TRUE,plotJunctions=NULL,nolegend=FALSE)
}
\arguments{
  \item{dfplot}{
the data.frame returned from the ViewAndScoreFull function
}
  \item{mismatchcolor}{
display mismatches as this color
}
\item{title}{
plot title
}
\item{covMaxDisplay}{
display at a maximum this number of read-pairs
}
\item{flipLeftandRight}{
flip the alignment views of the two loci so that the junction displayed
on the right will now be on the left and vice versa
}
\item{cropGaps}{
Only display bases that contain aligned reads
}
\item{plotJunctions}{
names of the two junctions plotted
}
\item{nolegend}{
do not show a legend
}}
\details{
see example
}
\value{
A list containing two items; (1) a ggplot object and (2) the names of
the plotted junctions.

}
\author{
Eitan Halper-Stromberg
}
\examples{
## this function is called by the plotSV method
}
