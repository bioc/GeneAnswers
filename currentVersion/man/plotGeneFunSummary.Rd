\name{plotGeneFunSummary}
\alias{plotGeneFunSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot the summarized gene annotation (ontologies) of a gene}
\description{
plot ontology graphs of the summarized gene annotation (ontologies), which is return by function \code{\link{geneFunSummarize}}
}
\usage{
plotGeneFunSummary(geneFunSummarizeResult, onto.graph, selGene = NULL, onto.graph.closure = NULL, allOntoID.direct = NULL, fdr.adjust = "none", showMiniSet = TRUE, highlightBest = TRUE, miniSetPvalue = 10^-5, ID2Name = NULL, savePrefix = "", geneSymbol = TRUE, saveImage = TRUE, selectionMethod = c("simple", "bestOnly", "all"), lib = "org.Hs.eg.db", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneFunSummarizeResult}{the output of \code{\link{geneFunSummarize}} function}
  \item{onto.graph}{an ontology graph in graphNEL class, which keeps the graph of the entire ontology DAG. Each node of it is an ontology ID}
  \item{selGene}{a subset of genes to plot ontology graph. If it is NULL, all plots of all genes in the "geneFunSummarizeResult" will be plotted.}
  \item{onto.graph.closure}{a flatten ontology graph in graphNEL class, whose edges represent offspring relationships.}
  \item{allOntoID.direct}{a vector of ontology IDs, which has direct gene mappings.}
  \item{fdr.adjust}{FDR adjustment methods used in p-value estimation.}
  \item{showMiniSet}{whether to high-light miniSet in the ontology plot}
  \item{highlightBest}{whether to high-light the most significant ontology term.}
  \item{miniSetPvalue}{threshold of miniSet p-value}
  \item{ID2Name}{a named vector of ontology IDs to ontology terms mapping}
  \item{savePrefix}{the filename prefix of the image files}
  \item{geneSymbol}{whether use gene symbol in the file name.}
  \item{saveImage}{whether save images as files.}
  \item{selectionMethod}{the method to select ontology nodes included in the plot graph}
  \item{lib}{the library used to get gene symbols of the Entrez Gene ID}
  \item{\dots}{other parameters used by function \code{\link{plotOntologyGraph}} }
}
\details{

}
\value{
Invisibly return a list of ontology subgraph (in graphGEL class), which correspond to the genes in the "geneFunSummarizeResult" (or "selGene" if provided as an input parameter).
}

\author{
Pan DU
}

\seealso{
See Also \code{\link{plotGraph}} and \code{\link{plotOntologyGraph}} 
}

\examples{

data(DO)
geneSummary <- geneFunSummarize('4267', gene2DO.map, DO.graph.closure.gene)[[1]]
plotGeneFunSummary(geneSummary['4267'], onto.graph=DO.graph.gene, onto.graph.closure=DO.graph.closure.gene, ID2Name=DO.terms, p.value.th=0.001, miniSetPvalue=10^-5, saveImage=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{hplot}
