\name{topCategoryGenes}
\alias{topCategoryGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Present top enrichment test information with genes }
\description{
  Function to present top enichmentInfo of given GeneAnswers instance with genes.
}
\usage{
topCategoryGenes(inputX, orderby = c("geneNum", "pvalue", "foldChange", "oddsRatio", "correctedPvalue"), top = 5, genesOrderBy = 1, decreasing = FALSE, topGenes = 5, file = FALSE, fileName = "topCategoryGenes.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputX}{ a given GeneAnswers instance }
  \item{orderby}{ type to sort enrichmentInfo slot }
  \item{top}{ integer to specify how many top rows to be presented }
  \item{genesOrderBy}{ integer or characters to specify gene ordered column }
  \item{decreasing}{ logic value to specify gene order is descending or not }
  \item{topGenes}{ integer to speicify how many top genes to be presented }
  \item{file}{ logic value to determine whether save to a file }
  \item{fileName}{ string to specify file name, default file name is topCategoryGenes.txt}
}
\details{
  orderby could be one of 'geneNum', 'pvalue', 'foldChange', 'oddsRatio' and 'correctedPvalue'.
  top could be an integer or 'ALL'. The top former specified categories will be printed on screen while only 30 categories will be displayed for 'ALL'. All categories can be saved in a specified file. topGenes is similar to top, but only top 5 genes will be displayed for 'ALL'.
  genesOrderBy could be an integer to specify column to be sorted. It can also be the column name.
}
\value{
  print necessary information on the screen and save into a specified file if request.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
# x is a GeneAnswers instance
\dontrun{topCategoryGenes(x, orderby='p')}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }

