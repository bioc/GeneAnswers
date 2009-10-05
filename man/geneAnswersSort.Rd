\name{geneAnswersSort}
\alias{geneAnswersSort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Sort enrichmentInfo of a GeneAnswers instance }
\description{
  a function to sort enrichmentInfo data frame in GeneAnswers objects.
}
\usage{
geneAnswersSort(x, sortBy = c("geneNum", "pvalue", "foldChange", "oddsRatio", "correctedPvalue"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a GeneAnswers instance }
  \item{sortBy}{ sorted type }
}
\details{
  sortBy could be one of "geneNum", "pvalue", "foldChange", "oddsRatio" and "correctedPvalue".
}
\value{
  return a new GeneAnswers instance with sorted by the specified type.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Gang Feng, Pan Du and Simon Lin}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
example(GeneAnswers)
xx <- geneAnswersSort(x, sortBy='correctedPvalue')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
