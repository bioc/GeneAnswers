\name{geneConceptNet}
\alias{geneConceptNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate Concept-gene network }
\description{
  Function to generate concept-gene network based on given list.
}
\usage{
geneConceptNet(inputList, lengthOfRoots=NULL, inputValue = NULL, centroidSize = "geneNum", output = c("fixed", "interactive", "none"), colorMap=NULL, bgColor='#ffffff', matchMode=c('absolute', 'relative'), zeroColorIndex=NULL, verbose=FALSE, symmetry=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputList}{ a character list to generate concept-gene network. Names of the list are concepts. }
  \item{lengthOfRoots}{ an integer, how many first elements could be root nodes. }
  \item{inputValue}{ NULL or a numeric vector to be used for color of nodes. }
  \item{centroidSize}{ 'geneNum' or a numeric vector to specify the size of concept nodes. }
  \item{output}{ type to specify output figure types }
  \item{colorMap}{ a R compatible color character vector, or NULL by embedded color scheme. }
  \item{bgColor}{ a R compatible color, default is '#ffffff' (white)}
  \item{matchMode}{ the mode of values matching colors, valid only if inputValue is not NULL, see details}
  \item{zeroColorIndex}{ index of color corresponding to zero, see details}
  \item{verbose}{ logic, determine whether show messages}
  \item{symmetry}{ logic, determine whether positive and negative values use the same color level.}
}
\details{
  The color of gene nodes could be specified by inputValue. Its length should be same as the total number of unique genes in inputList. There are two type of color matching methods. 'absolute' means, given zeroColorIndex that is color index in the colorMap for value 0, any value more than 0 will be matched to color between zeroColorIndex and the last one in colorMap based on the ratio of the value tothe maximum of the inputValue, while the value less than 0 will be matched to color between the first color in colorMap and zeroColorIndex, also based on the ratio of the value to the minimum of the inputValue. 'relative' means, set the first and last colors in colorMap to minimum and maximum of the inputValue, respectively, then any value between them will be mapped. If colorMap is set to NULL, the default color scheme will be applied. If the matching method is 'absolute', the color of 0 or the median of inputValue for 'relative' method, is set by bgColor, default value is '#ffffff' (white). The most positive value is represented as '#ff0000' (red), '#00ff00' (green) for the most negative value. 
  
  There are two types of output figures. "Fixed" means a network will be drawn on a regular R canvas, while "interactive" will generate a tck/tk canvas. Users can adjust nodes on it by mouse. "none" means no graphics output and return the attributes of vertices and edges.
}
\value{
  a concept-gene network is generated. A 3-element (1st one: igraph object; 2nd one: a dataframe for vertices attributes; 3rd one: a dataframe for edge attributes) list is returned when output is set to "none".
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\examples{
input <- list('ele01'=c('Aa', 'Bb'), 'ele02'=c('Bb', 'Cc', 'dd'))
\dontrun{geneConceptNet(input)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
