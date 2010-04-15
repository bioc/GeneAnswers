\name{topGOGenes}
\alias{topGOGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Present top GO enrichment test information with genes }
\description{
  Function to present top GO enichmentInfo of given GeneAnswers instance with genes. 
}
\usage{
topGOGenes(x, catTerm = TRUE, geneSymbol = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a given GeneAnswers instance with GO test }
  \item{catTerm}{ logic value to determine whether mapping GO IDs to GO terms }
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
##x is a GeneAnswers instance with GO test
\dontrun{topGOGenes(xxx, geneSymbol=F, catTerm=F, orderby='p')}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
