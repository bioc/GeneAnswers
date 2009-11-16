\name{buildNet}
\alias{buildNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ build and display a network for given IDs and interaction Matrix }
\description{
  A function to build and display a network for given IDs and interaction Matrix with specified filtered IDs.
}
\usage{
buildNet(graphIDs, idType=c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'GeneInteraction', 'Customized'), edgeM=NULL, layers=1, filterGraphIDs=NULL, filterLayer=0,
 					annLib=c('org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db', 'org.Dm.eg.db', 'customized'), output=c('interactive', 'fixed'), netMode=c('layer', 'connection'),
					vertexSize = NULL, edgeColor = NULL, colorMap=NULL, zeroColorIndex=NULL, matchMode=c('absolute', 'relative'),  label=TRUE, steric=FALSE, 
					directed=FALSE, direction=c('up', 'down', 'both'), showModeForNodes=c('nodes', 'filters'), verbose=TRUE, readable=TRUE, labelSize=1, labelColor='#666666',  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graphIDs}{ a character vector for given IDs }
  \item{idType}{ type of IDs, could be one of 'GO', 'GO.BP', 'GO.CC', 'GO.MF', 'GeneInteraction' and 'Customized'  }
  \item{edgeM}{ a 2-column Matrix representing a network  }
  \item{layers}{ an integer, specify how many layers will be retrieved. }
  \item{filterGraphIDs}{ a chacater vector for filtered IDs or a 2- or 3-column matrix for extra values.}
  \item{filterLayer}{ an integer, specify where filterGraphIDs are applied. }
  \item{annLib}{ type for annotation library, 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db', 'org.Dm.eg.db'and 'customized'. For 'customized', edgeM is necessary}
  \item{output}{ type to specify output figure types }
  \item{netMode}{ type to show network, see details}
  \item{vertexSize}{ an integer, the size of vertices in the network, default is NULL}
  \item{edgeColor}{ a R compatible color type, the color of edges in the network, default is NULL}
  \item{colorMap}{ a R compatible color character vector, or NULL by embedded color scheme. }
  \item{zeroColorIndex}{ index of color corresponding to zero, see details}
  \item{matchMode}{ the mode of values matching colors, valid only if inputValue is not NULL, see details}
  \item{label}{ logic, specify whether put labels for non-given nodes in the network. }
  \item{steric}{ logic, specify whether show 3D network. Because of igraph limited support for 3D, this function is just for fun!}
  \item{directed}{ logic, the network is a directed or not}
  \item{direction}{ search direction, it could be 'up', 'down' and 'both'. Valid for directed network only.}
  \item{showModeForNodes}{ type, the show mode for nodes on the network, only valid if filterGraphIDs is not NULL, see details}
  \item{verbose}{ logic, specify to show information or not.}
  \item{readable}{ logic, specify whether show IDs or Terms/Names for nodes}
  \item{labelSize}{ an integer, the size of label for nodes}
  \item{labelColor}{ an R compatible color, default is #666666}
  \item{\dots}{ other parameters used by 'getCategoryTerms' } 
}
\details{
  Currently, if idType is 'GO', 'GO.BP', 'GO.CC' or 'GO.MF', edgeM will be ignore.

  edgeM is a 2-column matrix. For directional connection, the direction is from column 1 elements to column 2 elements. For non-directional connection, each connection should be reversely presented twice, one is from column 1 element to column 2 element, while another is from column 2 element to column 1 element. In other words, non-directional connection is considered as two reverse directional connections.

  filterGraphIDs are applied only at the filterLayer and more outer layers. This means the nodes between the filterLayer layer and the most external layer belong to the filterGraphIDs. The nodes between given graphIDs and the (filterLayer-1) layer are or are not from filterGraphIDs, but those nodes not in filterGraphIDs should be able to be finally connected by given graphIDs and filterGraphIDs.

  There are two type of color matching methods. 'absolute' means, given zeroColorIndex that is color index in the colorMap for value 0, any value more than 0 will be matched to color between zeroColorIndex and the last one in colorMap based on the ratio of the value tothe maximum of the inputValue, while the value less than 0 will be matched to color between the first color in colorMap and zeroColorIndex, also based on the ratio of the value to the minimum of the inputValue.
  
  showModeForNodes stands for, if the filterGraphIDs is not NULL, some or all of filterGraphIDs could be nodes for given IDs multiple search. If it is set to 'nodes', it means only the values of nodes in the display network will be used to match color by matchMode. For 'filters', it means the values of all filter nodes will be used to match color. If values for color of nodes in the network are not large, while the maximum of color of filter nodes is large, it is recommended to set to 'nodes', or it is difficult to see difference for the nodes. For comparing two networks, for example, one is up-search and another is down-search for the same IDs, it is better to set to 'absolute' for easy comparisons.
 
  There are two types of output figures. "Fixed" means a network will be drawn on a regular R canvas, while "interactive" will generate a tck/tk canvas. Users can adjust nodes on it by mouse.

  If the filterGraphIDs is a ID vector. The filterGraphIDs nodes will be black, others will be white. If filterGraphIDs is a 2- or 3-column matrix, the 1st column is filter IDs and 2nd column is for color of nodes. If the 3rd column is available, it is for size of nodes.
  
  There are two types of netMode. 'layer' means size of nodes will be smaller and smaller for more and more external layers. And also color of edges change for different layers. 'connection' mode just distinguish direct or indirect connection. The size of the given IDs the largest. However, if filterGraphIDs is a 3-column matrix, the size of nodes will be determined by the 3rd column of filterGraphIDs.
  
  The graphIDs nodes are yellow circled solid dots. Color depends on colorMap and filterGraphIDs 2nd column. If no value available, all given graphIDs filterGraphIDs nodes are black, others are white.         
}
\value{
  invisibly return a list containing elements to represent a network. 
}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{getCategoryTerms}}, ~~~ }
\examples{
require(GeneAnswers)
example(GeneAnswers)
filterM <- cbind(rownames(getEnrichmentInfo(x)), -log2(getEnrichmentInfo(x)[,7]), getEnrichmentInfo(x)[,1])

\dontrun{buildNet(rownames(getEnrichmentInfo(x))[6:9], layers=5, filterGraphIDs=filterM, filterLayer=3, directed=TRUE, output='fixed')}
\dontrun{buildNet(rownames(getEnrichmentInfo(x))[242:244], layers=2, filterGraphIDs=filterM, filterLayer=1, directed=TRUE, output='fixed', netMode='connection', direction='down')}
\dontrun{buildNet(rownames(getEnrichmentInfo(x))[6:9], layers=3, filterGraphIDs=filterM[,1:2], filterLayer=3, directed=TRUE, output='fixed', netMode='connection')}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
