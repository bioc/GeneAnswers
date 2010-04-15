\name{geneAnswersReadable}
\alias{geneAnswersReadable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Make GeneAnswers Instance readable }
\description{
  a function to mapping category IDs and gene IDs to names and symbols.
}
\usage{
geneAnswersReadable(x, catTerm = TRUE, geneSymbol = TRUE, strict = FALSE, verbose=TRUE, missing=c('name', 'keep', 'remove'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a GeneAnswers instance containing category IDs and geneIDs }
  \item{catTerm}{ logic value to determine whether mapping category IDs to names }
  \item{geneSymbol}{ logic value to determine whether mapping gene IDs to symbols }
  \item{strict}{ logic value to determine whether interrupt conversion if NA is introduced. }
  \item{verbose}{ logical, show current stage or not}
  \item{missing}{ type of handling NA mapping. }
}
\details{
  Conversion could stop if NA is introduced and strict is set to TRUE.
  There are three types of parameters for variable 'missing'. 'name' means the NA mapping values are replaced by their names. 'keep' means all of NA values are kept. 'remove' means all of NA values are removed. 
}
\value{
  return a GeneAnswers instance with category names and/or gene symbols.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ \code{\link{getSymbols}}, \code{\link{getCategoryTerms}}}
\examples{
example(GeneAnswers)
xx <- geneAnswersReadable(x) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
