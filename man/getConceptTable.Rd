\name{getConceptTable}
\alias{getConceptTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate top concepts-genes table }
\description{
  Function to generate a top concepts-genes table based on a given GeneAnswers instance list.
}
\usage{
getConceptTable(gAList, topCat=10, items=c('both', 'geneNum', 'pvalue'),  sortBy = c('pvalue', 'geneNum', 'foldChange', 'oddsRatio', 'correctedPvalue'), catTerm=TRUE, strict=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gAList}{ a GeneAnswers instance list }
  \item{topCat}{ a numeric or string vector specified categories  } 
  \item{items}{ specify the contents in cells, see details }
  \item{sortBy}{ sorted type }
  \item{catTerm}{ a logic value to specify whether mapping category IDs to category names }
  \item{strict}{ logic value to stop conversion if NA is introduced. }
}
\details{
  A list containing two top concepts-genes tables is generated. The first table consists of gene amounts and enrichment test p values if 'items' is set to 'both'. 
Only gene amounts are kept if items is set to 'geneNum' or enrichment test p values if it is set to 'p values', while the second table contains enrichment test p values
}
\value{
  return a concepts-genes matrix list.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{\code{\link{geneAnswersBuilder}} }
\examples{
data(sampleGroupsData)
gAKEGGL <- lapply(sampleGroupsData, geneAnswersBuilder, 'org.Hs.eg.db', categoryType='KEGG', pvalueT=0.1, verbose=FALSE)
output<- getConceptTable(gAKEGGL)  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
