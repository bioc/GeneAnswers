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
\references{ ~put references to the literature/web site here ~ }
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data('DOLiteTerm')
getDOLiteTerms(c('DOLite:25', 'DOLite:142'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
