\name{categoryNet}
\alias{categoryNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Category Links }
\description{
  Function to plot a linkages of specified categories.
}
\usage{
categoryNet(catGenesList, centroidSize=NULL, output=c('fixed','interactive'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{catGenesList}{ a list of categories. }
  \item{centroidSize}{ a numeric vector to specify the size of concept nodes. If NULL, all of concept nodes are represented as the same size solid circles.}
  \item{output}{ type to specify output figure types. }
}
\details{
  catGenesList is a list of categories. Each element contains the genes in the corresponding category, respectively. And the names of the list are categories.
  If centroidSize is a numeric vector, its values are mapped to the categories in the catGenesList sequentially.
}
\value{
  A category linkage is generated.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
input <- list('cat1'=c(1,4,2,5), 'cat2'=c(3,5,8,9), 'cat3'=c(2,4,5,9), 'cat4'=c(1,5,3))
\dontrun{categoryNet(input)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
