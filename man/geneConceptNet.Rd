\name{geneConceptNet}
\alias{geneConceptNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate Concept-gene network }
\description{
  Function to generate concept-gene network based on given list.
}
\usage{
geneConceptNet(inputList, lengthOfRoots=NULL, inputValue = NULL, centroidSize = "geneNum", output = c("fixed", "interactive"), colorMap=NULL, matchMode=c('absolute', 'relative'), zeroColorIndex=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputList}{ a character list to generate concept-gene network. Names of the list are concepts. }
  \item{lengthOfRoots}{ an integer, how many first elements could be root nodes. }
  \item{inputValue}{ NULL or a numeric vector to be used for color of nodes. }
  \item{centroidSize}{ 'geneNum' or a numeric vector to specify the size of concept nodes. }
  \item{output}{ type to specify output figure types }
  \item{colorMap}{ a R compatible color character vector, or NULL by embedded color scheme. }
  \item{matchMode}{ the mode of values matching colors, valid only if inputValue is not NULL, see details}
  \item{zeroColorIndex}{ index of color corresponding to zero, see details}

}
\details{
  The color of gene nodes could be specified by inputValue. Its length should be same as the total number of unique genes in inputList. There are two type of color matching methods. 'absolute' means, given zeroColorIndex that is color index in the colorMap for value 0, any value more than 0 will be matched to color between zeroColorIndex and the last one in colorMap based on the ratio of the value tothe maximum of the inputValue, while the value less than 0 will be matched to color between the first color in colorMap and zeroColorIndex, also based on the ratio of the value to the minimum of the inputValue. 'relative' means, set the first and last colors in colorMap to minimum and maximum of the inputValue, respectively, then any value between them will be matched.
  
  There are two types of output figures. "Fixed" means a network will be drawn on a regular R canvas, while "interactive" will generate a tck/tk canvas. Users can adjust nodes on it by mouse.
}
\value{
  a concept-gene network is generated.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
input <- list('ele01'=c('Aa', 'Bb'), 'ele02'=c('Bb', 'Cc', 'dd'))
\dontrun{geneConceptNet(input)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
