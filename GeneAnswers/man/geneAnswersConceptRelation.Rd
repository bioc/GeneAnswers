\name{geneAnswersConceptRelation}
\alias{geneAnswersConceptRelation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Display a network related to given concepts for a GeneAnswers instance }
\description{
  A function to display a network related to given concepts of a GeneAnswer instance
}
\usage{
geneAnswersConceptRelation(x, showCats=c(1:5), conceptsIDs=NULL, directed=TRUE, direction=c('down', 'up', 'both'), catTerm=TRUE, catID=FALSE, nameLength='all', ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a GeneAnswers instance }
  \item{showCats}{ a numeric or string vector specified categories  }
  \item{conceptsIDs}{ a vector or a data frame or matrix containing possible relative concepts, see details}
  \item{directed}{ logic, the network is a directed or not}
  \item{direction}{ search direction, it could be 'up', 'down' and 'both'. Valid for directed network only.}
  \item{catTerm}{ a logic value to specify whether mapping category IDs to category names }
  \item{catID}{ a logic value to specify whether show category IDs when catTerm is set to TRUE }
  \item{nameLength}{ show how many first letters for long term names, 'all' for full name} 
  \item{\dots}{ other parameters used by 'getConnectedGraph' }
}
\details{
  conceptsIDs could be a character vector or a data frame or a matrix. As a character vector, it is a group of concept IDs or names depending on the given GeneAnswers instance, 
which are used to be a group of filters to draw a network relative to given concepts specified by showCats. When it is a data frame or matrix, it could be a 2- or 3-column data frame or matrix. The column 2 is always used to be represent nodes color, while the 3rd column is for size of nodes if available.
}
\value{
  return a invisible list representing the network. 
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ \code{\link{getConnectedGraph}} }
\examples{
require(GeneAnswers)
example(GeneAnswers)
\dontrun{geneAnswersConceptRelation(x, UP=FALSE, directed=TRUE, netMode='connection')}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
