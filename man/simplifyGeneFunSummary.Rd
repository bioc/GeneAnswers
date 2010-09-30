\name{simplifyGeneFunSummary}
\alias{simplifyGeneFunSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simplify the significant ontology terms to a miniSet annotation}
\description{Simplify the significant ontology terms to a mini-set, which includes the non-overlapping most significant terms and some other ontology terms, which have direct gene mapping but not included in the significant ontology terms.
}
\usage{
simplifyGeneFunSummary(geneFunSummarizeTestInfo, Onto.graph.closure, allOntoID.direct = NULL, fdr.adjust = "none", p.value.th = 10^-5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneFunSummarizeTestInfo}{the output of \code{\link{geneFunSummarize}} function}
  \item{Onto.graph.closure}{a graphNEL object, whose edges represent the link between a ontology term and its offspring ontology terms}
  \item{allOntoID.direct}{a vector of ontology IDs which has direct gene association. If not NULL, it will be used to filter the significant ontology terms}
  \item{fdr.adjust}{the FDR estimation methods used to estimate the FDR of function enrichment}
  \item{p.value.th}{the p-value threshold used to determine the significance of function enrichment}
}
\details{

}
\value{
The function return is a list with the same length of input "geneFunSummarizeTestInfo". Each element of the list is the simplified summarization of a gene. The simplified gene summarization is also a list with the following items
 \item{keptSigOntoID}{the significant ontology terms kept in the miniSet annotation}
 \item{keptEvidences}{the non significant ontology terms which has direct gene mapping but are not the offspring of the ontology terms in keptSigOntoID }
 \item{scores}{the annotation scores of all ontology terms in the miniSet annotation}
 The gene summarization object also includes the attributes of "pValueT".
}
\references{
Pan Du, Simon Lin, Gilbert Feng, Warren Kibbe, "GeneRIFcompendiate: Ranked gene annotations using collective GeneRIF associations and ontology terms", under review
}
\author{
Pan Du, Simon Lin, Gilbert Feng, Warren Kibbe
}

\seealso{
 See Also \code{\link{plotGeneFunSummary}} and \code{\link{geneFunSummarize}} 
}
\examples{
	data(DO)
	## test all ontology terms related with gene 643387
	geneSummary <- geneFunSummarize('643387', gene2DO.map, DO.graph.closure.gene, fdr.adjust = "none")
	## the p.values of all related ontology terms
	pValue <- sapply(geneSummary[[1]]$sigOntoInfo, function(x) x$pValue)
	pValue
	## simplify the annotation as a miniSet annotaiton
	geneSummary.sim <- simplifyGeneFunSummary(geneSummary, DO.graph.closure.gene, p.value.th=10^-5)
	geneSummary.sim
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
