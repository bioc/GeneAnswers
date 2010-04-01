\name{getConceptTable}
\alias{getConceptTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate top concepts-genes table }
\description{
  Function to generate a top concepts-genes table based on a given GeneAnswers instance list.
}
\usage{
getConceptTable(gAList, showCat=10, items=c('both', 'geneNum', 'pvalue'), catTerm=TRUE, strict=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gAList}{ a GeneAnswers instance list }
  \item{showCat}{ a numeric or string vector specified categories  } 
  \item{items}{ specify the contents in cells, see details }
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
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(sampleGroupsData)
gAKEGGL <- lapply(sampleGroupsData, geneAnswersBuilder, 'org.Hs.eg.db', categoryType='KEGG', pvalueT=0.05, FDR.correct=TRUE, verbose=FALSE)
output<- getConceptTable(gAKEGGL)  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
