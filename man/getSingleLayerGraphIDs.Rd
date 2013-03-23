\name{getSingleLayerGraphIDs}
\alias{getSingleLayerGraphIDs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ retrieve direct interacted nodes for given IDs and interaction Matrix }
\description{
  A function to retrieve direct interacted nodes for given IDs and interaction Matrix with specified filtered IDs.
}
\usage{
getSingleLayerGraphIDs(graphIDs, edgeM, remove=TRUE, filterGraphIDs=NULL, directed=FALSE, UP=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graphIDs}{ a character vector for given IDs }
  \item{edgeM}{ a 2-column Matrix representing connectionship  }
  \item{remove}{ logic, remove the non-connection graphIDs in the return values}
  \item{filterGraphIDs}{ a chacater vector for filtered IDs}
  \item{directed}{ logic, the network is a directed or not}
  \item{UP}{ logic, determine search Parents or Children. Only valid for directed relation.}
}
\details{
  edgeM is a 2-column matrix. For directional connection, the direction is from column 1 elements to column 2 elements. For 
non-directional connection, each connection should be reversely presented twice, one is from column 1 element to column 2 element, while
another is from column 2 element to column 1 element. In other words, non-directional connection is considered as two reverse
directional connections.
  filterGraphIDs is used to only keep nodes in filterGraphIDs.  
}
\value{
  return a list representing a network.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\examples{
m <- matrix(c('1','4', '2', '6', '1', '5', '3', '7', '5', '2'), ncol=2, byrow=TRUE)
m
getSingleLayerGraphIDs(c('1','2','3'), m)

# if the connection is not directional, the connection between '5' and '2' will be missed without changing m.
m <- rbind(m, c('2', '5'))
getSingleLayerGraphIDs(c('1','2','3'), m)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
