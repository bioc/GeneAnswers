\name{topDOLiteGenes}
\alias{topDOLiteGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Present top DOLite enrichment test information with genes}
\description{
  Function to present top DOLite enichmentInfo of given GeneAnswers instance with genes. 
}
\usage{
topDOLiteGenes(x, catTerm = TRUE, geneSymbol = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a given GeneAnswers instance with DOLite test }
  \item{catTerm}{ logic value to determine whether mapping DOLite IDs to DOLite terms }
  \item{geneSymbol}{ logic value to determine whether mapping gene Entrez IDs to gene symbols}
  \item{\dots}{ other parameters to transfer to topCategoryGenes }
}
\details{
  See function topCategoryGenes help for details
}
\value{
  print necessary information on the screen and save into a specified file if request.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{\code{\link{topCategoryGenes}}}
\examples{
##x is a GeneAnswers instance with DOLite test
\dontrun{topDOLiteGenes(x, geneSymbol=TRUE, orderby='pvalue', top=10, topGenes='ALL', genesOrderBy='pValue', file=TRUE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }

