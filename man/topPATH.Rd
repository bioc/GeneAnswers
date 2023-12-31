\name{topPATH}
\alias{topPATH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Present top KEGG enrichment test information }
\description{
  Function to present top KEGG enichmentInfo of given GeneAnswers instance.
}
\usage{
topPATH(x, catTerm = TRUE, keepID = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a given GeneAnswers instance containing KEGG information }
  \item{catTerm}{ logic value to determine whether mapping to DOLite terms or not }
  \item{keepID}{ logic value to determine whether showing IDs or not}
  \item{\dots}{ other parameters to transfer to topCategory }
}
\value{
  print necessary information on the screen and save into a specified file if request.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{\code{\link{topCategory}}}
\examples{
# x is a GeneAnswers instance with DOLite test
\dontrun{topPATH(x, top=10)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }

