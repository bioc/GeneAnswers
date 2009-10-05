\name{geneConceptNet}
\alias{geneConceptNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate Concept-gene network }
\description{
  Function to generate concept-gene network based on given list.
}
\usage{
geneConceptNet(inputList, inputValue = NULL, centroidSize = "geneNum", output = c("fixed", "interactive"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputList}{ a character list to generate concept-gene network. Names of the list are concepts. }
  \item{inputValue}{ NULL or a numeric vector to be used for color of nodes. }
  \item{centroidSize}{ 'geneNum' or a numeric vector to specify the size of concept nodes. }
  \item{output}{ type to specify output figure types }
}
\details{
  The color of gene nodes could be specified by inputValue. Its length should be same as the total number of unique genes in inputList. Genes with positive values will be represented by red solid circles, while green nodes stand for negative values gene nodes.
  There are two types of output figures. "Fixed" means a network will be drawn on a regular R canvas, while "interactive" will generate a tck/tk canvas. Users can adjust nodes on it by mouse.  
}
\value{
  a concept-gene network is generated.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
input <- list('ele01'=c('Aa', 'Bb'), 'ele02'=c('Bb', 'Cc', 'dd'))
\dontrun{geneConceptNet(input)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
