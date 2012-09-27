\name{saveGeneFunSummary}
\alias{saveGeneFunSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{save the summarized gene function as a tab-separated text file}
\description{save the summarized gene function as a tab-separated text file}
\usage{
saveGeneFunSummary(geneFunSummarizeResult, simplifyInfo = NULL, addFDR = FALSE, species = "human", ID2Name = NULL, fileName = "geneSummarization.xls")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneFunSummarizeResult}{the output of \code{\link{geneFunSummarize}} function}
  \item{simplifyInfo}{the output of \code{\link{simplifyGeneFunSummary}} function}
  \item{addFDR}{whether add FDR information or not}
  \item{species}{the species of the related ontology under test, which is used to convert the Entrez Gene ID as gene symbols}
  \item{ID2Name}{a named vector of ontology IDs to ontology terms mapping}
  \item{fileName}{the file name to keep the summarized gene information}
}
\details{

}
\value{
Invisibly return a data.frame of summarized gene function, which includes "geneID", "geneSymbol", "bestOntology", "enrichedOntology" and "allEvidence". The corresponding annotation scores are kept in the parenthesis behind the ontology IDs.
}
\references{

}
\author{
Pan DU
}

\seealso{
See also \code{\link{geneFunSummarize}}
}
\examples{
	data(DO)
	geneSummary <- geneFunSummarize(c("5037","9314"), gene2DO.map, DO.graph.closure.gene)
	summaryInfo <- saveGeneFunSummary(geneSummary, fileName='geneSummarization.xls')
	summaryInfo
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
