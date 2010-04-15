\name{topPATHGenes}
\alias{topPATHGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Present top KEGG enrichment test information with genes }
\description{
 Function to present top KEGG enichmentInfo of given GeneAnswers instance with genes.
}
\usage{
topPATHGenes(x, catTerm = TRUE, geneSymbol = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a given GeneAnswers instance with KEGG test }
  \item{catTerm}{ logic value to determine whether mapping KEGG IDs to KEGG terms }
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
\seealso{ ~~objects to See Also as \code{\link{topCategoryGenes}}, ~~~ }
\examples{
##x is a GeneAnswers instance with KEGG test
\dontrun{topPATHGenes(x, geneSymbol=TRUE, orderby='genenum', top=6, topGenes=8, genesOrderBy='foldChange')} 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
