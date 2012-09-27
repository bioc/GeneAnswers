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
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin}
\seealso{\code{\link{GeneAnswers-class}}}
\examples{
example(GeneAnswers)
xx <- geneAnswersSort(x, sortBy='correctedPvalue')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
