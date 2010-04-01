\name{topDOLite}
\alias{topDOLite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Present top DOLite enrichment test information }
\description{
  Function to present top DOLite enichmentInfo of given GeneAnswers instance.
}
\usage{
topDOLite(x, catTerm = TRUE, keepID=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a given GeneAnswers instance containing DOLite information }
  \item{catTerm}{ logic value to determine whether mapping to DOLite terms or not }
  \item{keepID}{ logic value to determine whether showing IDs or not}
  \item{\dots}{ other parameters to transfer to topCategory }
}
\value{
  print necessary information on the screen and save into a specified file if request.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
# x is a GeneAnswers instance with DOLite test
\dontrun{topDOLite(x, top=10)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }

