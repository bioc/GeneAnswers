\name{geneAnnotationHeatmap}
\alias{geneAnnotationHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Make a concept-gene cross tabulation }
\description{
   Function to make a concept-gene cross tabulation
}
\usage{
geneAnnotationHeatmap(annotationList, dataMatrix = NULL, addGeneLabel = TRUE, colorMap = c("#000000", "#FFFFFF"), sortBy = "both", standardize.data = TRUE, colorMap.data = "GBR", sortBy.data = "row", mar = c(1, 1, 8, 6), cex.axis = c(0.8, 0.8), mapType = c("table", "heatmap"), displayAll=FALSE, symmetry=FALSE, colorBar=FALSE, colorBarLabel=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annotationList}{ a list of annotation to gene mapping. }
  \item{dataMatrix}{ a 2-dimensional numeric matrix. If it is provided, it will be plot side by side with the annotation heatmap. }
  \item{addGeneLabel}{ logic, indicate whether add gene labels }
  \item{colorMap}{ vector to specify color map of the two-color annotation heatmap }
  \item{sortBy}{ string to specify whether to sort the annotation matrix by row, column, both row and column or none of them   }
  \item{standardize.data}{ logic, specify whether to standardize the dataMatrix by row~~ }
  \item{colorMap.data}{ string to specify color map of the dataMatrix heatmap }
  \item{sortBy.data}{ string to specify whether to sort the dataMatrix by row, column, both row and column or none of them }
  \item{mar}{ integer vector to speicify margin of the plot }
  \item{cex.axis}{ integer vector to specify the character size of row and column labels }
  \item{mapType}{string to specify concept-gene map type}
  \item{displayAll}{ logic, specify to show all of gene expression profile or remove redundant entries.}
  \item{symmetry}{ logic, indicate the values corresponding to two extreme colors are same if TURE. }
  \item{colorBar}{ logic, show colorbar or not}
  \item{colorBarLabel}{ character vector to show color bar label.}
}
\details{
  This function basically generates two maps in one canvas. Left side is a heatmap based on given expression matrix. Right side is a concept-gene map, which could be represented as two-color heatmap or table, depends on parameter "mapType". 
}
\value{
  The function will generate a map without return value.
}
\author{ Pan Du, Gang Feng and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
a <- list(group1 = c('a','b','c','d','f'), group2= c('b','d','e','a','g','h'))
b <- matrix(rnorm(48), nrow=8,ncol=6)
rownames(b) <- tolower(LETTERS[1:8])
colnames(b) <- c('ctrl1', 'ctrl2', 'ctrl3', 'treat1', 'treat2', 'treat3')
\dontrun{geneAnnotationHeatmap(a,dataMatrix=b)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

