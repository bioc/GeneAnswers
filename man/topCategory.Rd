\name{topCategory}
\alias{topCategory}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Present top enrichment test information }
\description{
  Function to present top enichmentInfo of given GeneAnswers instance.
}
\usage{
topCategory(inputX, orderby = c("geneNum", "pvalue", "foldChange", "oddsRatio", "correctedPvalue"), top = 5, file = FALSE, fileName = "topCategory.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputX}{ a given GeneAnswers instance }
  \item{orderby}{ type to sort enrichmentInfo slot }
  \item{top}{ integer to specify how many top rows to be presented }
  \item{file}{ logic value to determine whether save to a file }
  \item{fileName}{ string to specify file name, default file name is topCategory.txt}
}
\details{
  orderby could be one of 'geneNum', 'pvalue', 'foldChange', 'oddsRatio' and 'correctedPvalue'.
  top could be an integer or 'ALL'. The top former specified categories will be printed on screen while only 30 categories will be displayed for 'ALL'. All categories can be saved in a specified file.
}
\value{
  print necessary information on the screen and save into a specified file if request.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
# x is a GeneAnswers instance
\dontrun{topCategory(x, orderby='pvalue')}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
