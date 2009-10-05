\name{chartPlots}
\alias{chartPlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Pie Chart and Bar Plots }
\description{
  Make pie chart and bar plot based on the given data frame.
}
\usage{
chartPlots(x, chartType = c("pieChart", "barPlot", "all"), specifiedCols = c("genes in Category"), top = 5, newWindow=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a data frame to be used for pie chart and box plot }
  \item{chartType}{ plot type, "pieChart", "barPlot" or both could be specified. }
  \item{specifiedCols}{ the column will be used to be represented. }
  \item{top}{ number to specify how many first categories will be drawn. }
  \item{newWindow}{ logic, determine whether draw on a new canvas. }
  \item{\dots}{ additional arguments passed to piechart or barplot.}
}
\details{
  chartType could be pie chart, bar plot or both (parameter is "all").
  specifiedCols is the column that will be used to plot. It could be column name or number.
  If chartType is set to 'all', the barplot will be drawn on a new canvas whatever newWindow is set to TRUE or FALSE.
}
\value{
  A pie chart and/or barplot are generated depends on specification.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
x <- matrix(c(6,9,3,30,13,2,15,20), nrow = 4, ncol=2, byrow=FALSE,
               dimnames = list(c("group1", "group2", "group3", "group4"),
                               c("value1", "value2")))
chartPlots(x, chartType='all', specifiedCol = "value2", top = 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
