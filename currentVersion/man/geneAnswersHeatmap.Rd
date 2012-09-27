\name{geneAnswersHeatmap}
\alias{geneAnswersHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate Concept-Gene Tabulates }
\description{
  A function to generate specified Concept-Gene Tabulates
}
\usage{
geneAnswersHeatmap(x, showCats = c(1:5), catTerm = FALSE, geneSymbol = FALSE, catID=FALSE, nameLength='all', showAllGenes=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an instance of GeneAnswers objects }
  \item{showCats}{ a numeric or string vector specified categories  }
  \item{catTerm}{ a logic value to specify whether mapping category IDs to category names }
  \item{geneSymbol}{ a logic value to specify whether mapping gene IDs to gene symbols }
  \item{catID}{ a logic value to specify whether show category IDs when catTerm is set to TRUE }
  \item{nameLength}{ show how many first letters for long term names, 'all' for full name}
  \item{showAllGenes}{logic, show all genes in the heatmap or not} 
  \item{\dots}{ other parameters used by geneAnnotationHeatmap }
}
\details{
  This function generates concept-gene tabulates for an input GeneAnswers instance. The concept-gene tabulates contain two maps. Left side is a heatmap based on given expression matrix. Right side is a concept-gene map, which could be represented as two-color heatmap or table. 
}
\value{
  The function will generate a map without return value.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ \code{\link{geneAnnotationHeatmap}}}
\examples{
example(GeneAnswers)
\dontrun{geneAnswersHeatmap(x, catTerm=TRUE, geneSymbol=TRUE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
