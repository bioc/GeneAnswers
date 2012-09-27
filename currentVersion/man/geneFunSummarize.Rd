\name{geneFunSummarize}
\alias{geneFunSummarize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarize gene functions (annotations) based collective annotation evidences associated with ontology terms}
\description{Summarize gene functions (annotations) based collective annotation evidences associated with ontology terms}
\usage{
geneFunSummarize(genes, gene2Onto, Onto2offspring, rmOntoID = c("DOID:4", "DOID:63"), p.value.th = 0.01, fdr.adjust = "fdr", minNumTh = 2, includeTestOnto = TRUE, directMapConstraint = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genes}{a vector of Entrez Gene IDs to do gene function summarization}
  \item{gene2Onto}{a list gene 2 Ontology mapping (Ontology IDs, duplicated IDs are allowed, which is equivalent to multiple  evidences with the same function)}
  \item{Onto2offspring}{a list or graphNEL object shows the relations of a ontology id to all its offsprings}
  \item{rmOntoID}{some ontology ids (like root id or too general ontology ids) can be pre-removed in the estimation}
  \item{p.value.th}{the p-value threshold used to determine the significance of function enrichment}
  \item{fdr.adjust}{the FDR estimation methods used to estimate the FDR of function enrichment}
  \item{minNumTh}{the minimum number of evidences required to claim as significant enriched ontology term}
  \item{includeTestOnto}{ whether include the direct evidence of the testing ontology term itself}
  \item{directMapConstraint}{ whether only consider the ontology ids with direct evidence mappings }
}
\details{

}
\value{
The function return is a list with the same length of input "genes". Each element of the list is the summarization of a gene. The gene summarization is also a list with the following items
 \item{allEvidence}{all evidences (ontology terms) related with the testing gene}
 \item{sigOntoInfo}{the significant ontology terms associated with the testing gene}
 \item{bestOntoInfo}{the information of the most significant ontology term associated with the testing gene}
 The gene summarization object also includes the attributes of the input parameter settings.
}
\references{
Pan Du, Simon Lin, Gilbert Feng, Warren Kibbe, "GeneRIFcompendiate: Ranked gene annotations using collective GeneRIF associations and ontology terms", under review
}
\author{
Pan Du, Simon Lin, Gilbert Feng, Warren Kibbe
}

\seealso{
 See Also \code{\link{plotGeneFunSummary}} and \code{\link{simplifyGeneFunSummary}} 
}
\examples{
	data(DO)
	## test all ontology terms related with gene PEBP1 (Entrez Gene ID: 5037)
	geneSummary <- geneFunSummarize('5037', gene2DO.map, DO.graph.closure.gene, p.value.th=0.01, fdr.adjust='none')
	## the p.values of all related ontology terms
	pValue <- sapply(geneSummary[[1]]$sigOntoInfo, function(x) x$pValue)
	pValue
	## plot the relations of the summarized gene annotation
	plotGeneFunSummary(geneSummary, onto.graph=DO.graph.gene, onto.graph.closure=DO.graph.closure.gene, ID2Name=DO.terms, p.value.th=0.0001, miniSetPvalue=10^-5, saveImage=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
