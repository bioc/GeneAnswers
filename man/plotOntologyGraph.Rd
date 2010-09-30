\name{plotOntologyGraph}
\alias{plotOntologyGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot the ontology graph}
\description{plot the ontology graph in graphNEL class, the color of ontology is based on the enrichment p.value}
\usage{
plotOntologyGraph(onto.pValue, relatedOntoID, onto.graph, bestOntoID = NULL, onto.graph.closure = NULL, rootID = "DOID:4", ID2Name = NULL, p.value.th = 0.01, fillColor = "white", colorLevel = seq(2, 20, by = 2), relative.color = TRUE, fontsize = 15, colorMap = colorRampPalette(c("white", "red"))(length(colorLevel) + 1), selectionMethod = c("simple", "bestOnly", "all"), omitNode = TRUE, saveImageName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{onto.pValue}{a vector ontology p.values, whose names are corresponding ontology IDs}
  \item{relatedOntoID}{a vector of related ontology IDs, which will be highlighted.}
  \item{onto.graph}{an ontology graph in graphNEL class, which keeps the graph of the entire ontology DAG. Each node of it is an ontology ID}
  \item{bestOntoID}{the ontology ID which has most significant p-value}
  \item{onto.graph.closure}{a flatten ontology graph in graphNEL class, whose edges represent offspring relationships.}
  \item{rootID}{the root ontology ID in the ontology graph}
  \item{ID2Name}{a named vector of ontology IDs to ontology terms mapping}
  \item{p.value.th}{enrichment p.value threshold to claim as significant}
  \item{fillColor}{the default fill color of the ontology nodes}
  \item{colorLevel}{the level of colors to show the significant enriched ontology terms.}
  \item{relative.color}{determine whether  the most significant enriched node has the deepest color.}
  \item{fontsize}{the font size of ontology terms}
  \item{colorMap}{the color map to show the significant enriched ontology terms}
  \item{selectionMethod}{the method to select ontology nodes included in the plot graph}
  \item{omitNode}{whether to omit the intermediate insignificant ontology terms between the significant terms and terms with direct gene mapping}
  \item{saveImageName}{the name of image to be save. If it is NULL, then the image will not be saved.}
}
\details{

}
\value{
Invisibly return a list, which include following elements:
	\item{graph}{a graphNEL object includes the ontology subgraph to plot}
	\item{layoutAttrs}{layoutAttrs used by the Rgraphviz function \code{\link[Rgraphviz]{layoutGraph}}}
	\item{layoutNodeAttrs}{layoutNodeAttrs used by the Rgraphviz function \code{\link[Rgraphviz]{layoutGraph}}}
	\item{layoutEdgeAttrs}{layoutEdgeAttrs used by the Rgraphviz function \code{\link[Rgraphviz]{layoutGraph}}}
	\item{nodeRenderAttrs}{nodeRenderAttrs used by the graph function \code{\link[graph]{nodeRenderInfo}}}
	\item{edgeRenderAttrs}{edgeRenderAttrs used by the graph function \code{\link[graph]{edgeRenderInfo}}}
These invisibly return information can be used as the input of function \code{\link{plotGraph}}.
}

\author{
Pan DU
}

\seealso{
See Also \code{\link{plotGraph}} and \code{\link{plotGeneFunSummary}} 
}
\examples{
	data(DO)
	geneSummary <- geneFunSummarize('4267', gene2DO.map, DO.graph.closure.gene)[[1]]
	pValue <- sapply(geneSummary$sigOntoInfo, function(x) x$pValue)
	plotOntologyGraph(pValue, geneSummary$allEvidence, DO.graph.gene, onto.graph.closure=DO.graph.closure.gene, ID2Name=DO.terms, p.value.th=0.001)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{hplot}

