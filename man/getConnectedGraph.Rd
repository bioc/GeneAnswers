\name{getConnectedGraph}
\alias{getConnectedGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ build and display a network for given IDs }
\description{
  A function to build and display a network by different show types for given IDs and interaction Matrix with specified filtered IDs.
}
\usage{
getConnectedGraph(graphIDs, idType=c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'GeneInteraction', 'Customized'), edgeM=NULL, limitedLayers=FALSE, layers=6, treeMergeFilter=FALSE,
						searchAll=FALSE, showAllNodes=FALSE, directed=FALSE, direction=c('up', 'down', 'both'), filterGraphIDs=NULL, filterLayer=1, verbose=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graphIDs}{ a character vector for given IDs }
  \item{idType}{ type of IDs, could be one of 'GO', 'GO.BP', 'GO.CC', 'GO.MF', 'GeneInteraction' and 'Customized'  }
  \item{edgeM}{ a 2-column Matrix representing a network  }
  \item{limitedLayers}{ logic, user specified layers to stop search}
  \item{layers}{ an integer, specify how many layers will be retrieved. }
  \item{treeMergeFilter}{ logic, determine whether apply filterGraphIDs during searching a merged tree, see details}
  \item{searchAll}{ logic, determine whether search all nodes}
  \item{showAllNodes}{ logic, determine whether show all nodes based on searching result}
  \item{directed}{ logic, the network is a directed or not}
  \item{direction}{ search direction, it could be 'up', 'down' and 'both'. Valid for directed network only.}
  \item{filterGraphIDs}{ a chacater vector for filtered IDs or a 2- or 3-column matrix for extra values.}
  \item{filterLayer}{ an integer, specify where filterGraphIDs are applied. }
  \item{verbose}{ logic, specify to show information or not.}
  \item{\dots}{ other parameters used by 'buildNet' } 
}
\details{
  Currently, if idType is 'GO', 'GO.BP', 'GO.CC' or 'GO.MF', edgeM will be ignore.

  edgeM is a 2-column matrix. For directional connection, the direction is from column 1 elements to column 2 elements. For non-directional connection, each connection should be reversely presented twice, one is from column 1 element to column 2 element, while another is from column 2 element to column 1 element. In other words, non-directional connection is considered as two reverse directional connections.

  filterGraphIDs are applied only at the filterLayer and more outer layers. This means the nodes between the filterLayer layer and the most external layer belong to the filterGraphIDs. The nodes between given graphIDs and the (filterLayer-1) layer are or are not from filterGraphIDs, but those nodes not in filterGraphIDs should be able to be finally connected by given graphIDs and filterGraphIDs.

  The function at first searches a merged tree based on given IDs. During searching, filterGraphIDs could be applied if 'treeMergeFilter' is set to TRUE. If a merged tree is found, searching process stops unless 'searchAll' is set to TRUE. However, 'limitedLayers' is set to TRUE, searching process also stops when searching layers reach 'layers'. Only all filterGraphIDs specified nodes as well as given nodes will be displayed if 'showAllNodes' is set to FALSE, or all connected nodes will be displayed.

  See buildnet for network layout.
}
\value{
  invisibly return a list containing elements to represent a network. 
}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{buildNet}}, ~~~ }
\examples{
require(GeneAnswers)
example(GeneAnswers)
filterM <- cbind(rownames(getEnrichmentInfo(x)), -log2(getEnrichmentInfo(x)[,7]), getEnrichmentInfo(x)[,1])
\dontrun{getConnectedGraph(rownames(getEnrichmentInfo(x))[c(1:5)],  filterGraphIDs=filterM,  output='fixed', directed=TRUE, direction='up', showAllNodes=T, searchAll=T)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
