\name{geneAnswersChartPlots}
\alias{geneAnswersChartPlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Make pie chart and bar plot }
\description{
  Make pie chart and bar plot for given GeneAnswers instance
}
\usage{
geneAnswersChartPlots(x, chartType=c('pieChart', 'barPlot', 'all'), sortBy = c('geneNum', 'pvalue', 'foldChange', 'oddsRatio', 'correctedPvalue'), newWindow=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a GeneAnswers instance }
  \item{chartType}{ plot type, "pieChart", "barPlot" or both could be specified. }
  \item{sortBy}{ the column will be used to be represented. }
  \item{newWindow}{ logic, determine whether draw on a new canvas. }
  \item{\dots}{ additional arguments passed to piechart or barplot.}
}
\details{
  chartType could be pie chart, bar plot or both (parameter is "all").
  specifiedCols is the column of enrichmentInfo that will be used to plot. It could be one of 'genes in Category', 'p value' or 'fdr p value'.
  If chartType is set to 'all', the barplot will be drawn on a new canvas whatever newWindow is set to TRUE or FALSE. 
}
\value{
  A pie chart and/or barplot are generated depends on specification.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
example(GeneAnswers)
\dontrun{geneAnswersChartPlots(x)} 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
