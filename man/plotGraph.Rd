\name{plotGraph}
\alias{plotGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot and render a graphNEL object}
\description{plot and render a graphNEL object}
\usage{
plotGraph(graph, layoutAttrs = NULL, layoutNodeAttrs = NULL, layoutEdgeAttrs = NULL, nodeRenderAttrs = NULL, edgeRenderAttrs = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{graph}{a graphNEL object includes the ontology subgraph to plot}
 \item{layoutAttrs}{layoutAttrs used by the Rgraphviz function \code{\link[Rgraphviz]{layoutGraph}}}
 \item{layoutNodeAttrs}{layoutNodeAttrs used by the Rgraphviz function \code{\link[Rgraphviz]{layoutGraph}}}
 \item{layoutEdgeAttrs}{layoutEdgeAttrs used by the Rgraphviz function \code{\link[Rgraphviz]{layoutGraph}}}
 \item{nodeRenderAttrs}{nodeRenderAttrs used by the graph function \code{\link[graph]{nodeRenderInfo}}}
 \item{edgeRenderAttrs}{edgeRenderAttrs used by the graph function \code{\link[graph]{edgeRenderInfo}}}
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
See Also \code{\link{plotOntologyGraph}}
}
\examples{
	library(graph)
	data(DO)
	geneSummary <- geneFunSummarize('4267', gene2DO.map, DO.graph.closure.gene)[[1]]
	pValue <- sapply(geneSummary$sigOntoInfo, function(x) x$pValue)
	tmp <- plotOntologyGraph(pValue, geneSummary$allEvidence, DO.graph.gene, onto.graph.closure=DO.graph.closure.gene, ID2Name=DO.terms, p.value.th=0.001)
	# modify the graph
	tmpGraph <- tmp$graph
	nn <- nodes(tmpGraph)
	# remove a node in the graph
	tmpGraph <- subGraph(nn[nn != 'DOID:3347'], tmpGraph)
	tmpGraph <- addEdge("DOID:1115", "DOID:3376", tmpGraph, 2)
	edgeRenderAttrs <- tmp$edgeRenderAttrs
	edgeRenderAttrs$lty <- c(edgeRenderAttrs$lty, "DOID:1115~DOID:3376"="dashed")
	plotGraph(tmpGraph, layoutAttrs=tmp$layoutAttrs, nodeRenderAttrs=tmp$nodeRenderAttrs, edgeRenderAttrs=edgeRenderAttrs, layoutNodeAttrs=tmp$layoutNodeAttrs, layoutEdgeAttrs=tmp$layoutEdgeAttrs) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{hplot}
