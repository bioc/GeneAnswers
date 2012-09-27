\name{geneAnswersConcepts}
\alias{geneAnswersConcepts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Concept-Gene Networking Plotting }
\description{
  A function to generate a concept-gene network by given gene information
}
\usage{
geneAnswersConcepts(x, centroidSize=c('geneNum', 'pvalue', 'foldChange', 'oddsRatio', 'correctedPvalue'), output=c('fixed','interactive'), showCats=c(1:5), catTerm=FALSE, catID=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a GeneAnswers instance. }
  \item{centroidSize}{ type to represent the size of concepts. }
  \item{output}{ output type of final output. }
  \item{showCats}{ a numeric or string vector specified categories  }
  \item{catTerm}{ a logic value to specify whether mapping category IDs to category names }
  \item{catID}{ a logic value to specify whether show category IDs when catTerm is set to TRUE }
}
\details{
  centroidSize could be one of "geneNum", "pvalue", "foldChange", "oddsRatio", "correctedPvalue". Each one defines to which the size of cencept dot is proportional
  geneNum: number of genes connecting to the concept
  pvalue: p value of enrichment test
  foldChange: fold of gene overrepresent in concepts
  oddsRatio: odds ratio of enrichment test	
  correctedPvalue: adjusted p value of enrichment test
  output defines whether the final figure is interactive or not. Interactive figure calls igraph package to generate a tck/tk canvas. Fixed figure is a non-interactive png figure.
}
\value{
  One category-linkage figure is generated. It could be a R figure or tcltk figure depends on how the user set parameter output. 
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ \code{\link{categoryNet}} }
\examples{
example(GeneAnswers)
\dontrun{geneAnswersConcepts(x, centroidSize='pvalue', output='interactive')}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
