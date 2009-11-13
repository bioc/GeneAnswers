\name{getMultiLayerGraphIDs}
\alias{getMultiLayerGraphIDs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ retrieve multilayer interacted nodes for given IDs and interaction Matrix }
\description{
  A function to retrieve multilayer interacted nodes for given IDs and interaction Matrix with specified filtered IDs.
}
\usage{
getMultiLayerGraphIDs(graphIDs, idType=c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'GeneInteraction', 'Customized'), edgeM=NULL, layers=1, filterGraphIDs=NULL, filterLayer=0, UP=TRUE, directed=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graphIDs}{ a character vector for given IDs }
  \item{idType}{ type of IDs, could be one of 'GO', 'GO.BP', 'GO.CC', 'GO.MF', 'GeneInteraction' and 'Customized'  }
  \item{edgeM}{ a 2-column Matrix representing a network  }
  \item{layers}{ an integer, specify how many layers will be retrieved. }
  \item{filterGraphIDs}{ a chacater vector for filtered IDs}
  \item{filterLayer}{ an integer, specify where filterGraphIDs are applied. }
  \item{UP}{ logic, determine search Parents or Children. Only valid for directed relation.}
  \item{directed}{ logic, the network is a directed or not}
}
\details{
  Currently, if idType is 'GO', 'GO.BP', 'GO.CC' or 'GO.MF', edgeM will be ignore.
  edgeM is a 2-column matrix. For directional connection, the direction is from column 1 elements to column 2 elements. For 
non-directional connection, each connection should be reversely presented twice, one is from column 1 element to column 2 element, while
another is from column 2 element to column 1 element. In other words, non-directional connection is considered as two reverse
directional connections.
  filterGraphIDs are applied only at the filterLayer and more outer layers. This means the nodes between the filterLayer layer and the most external layer belong to the filterGraphIDs. The nodes between given graphIDs and the (filterLayer-1) layer are or are not from filterGraphIDs, but those nodes not in filterGraphIDs should be able to be finally connected by given graphIDs and filterGraphIDs.       
}
\value{
  return a list containing elements to represent a network. The first element is a logic value, TRUE means no more connection between the most external layer nodes and other nodes. The second element is a list of layer-length. If the 1st element is FASLE, the length of 2nd element should be (layers + 1). And starting from the 3rd elements, the remaining elements construct a network.
}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{getSingleLayerGraphIDs}}, ~~~ }
\examples{
require(GeneAnswers)
example(GeneAnswers)
getMultiLayerGraphIDs(rownames(getEnrichmentInfo(x))[5:6], UP=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
