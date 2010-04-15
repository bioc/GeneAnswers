\name{getDOLiteTerms}
\alias{getDOLiteTerms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get DOLite Terms of Given DOLite IDs }
\description{
  function to convert DOLite IDs to DOLite Terms
}
\usage{
getDOLiteTerms(DOLiteIDs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DOLiteIDs}{ a character vector containing DOLite IDs }
}
\value{
  return a DOLite term vector based on given DOLite IDs.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ \code{\link{getCategoryTerms}}}
\examples{
data('DOLiteTerm')
getDOLiteTerms(c('DOLite:25', 'DOLite:142'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
