\name{getcaBIOPATHList}
\alias{getcaBIOPATHList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Retrieve caBIO path categories containing given genes }
\description{
  Function to retrieve caBIO pathway IDss containing the given genes.
}
\usage{
getcaBIOPATHList(lls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lls}{ an Entrez gene IDs vector }
}
\details{
  The given gene IDs should be Entrez gene IDs. And the return list also only contains Entrez gene IDs besides caBIO pathway IDs.  
}
\value{
  return an Entrez genes ID list, names of the list are caBIO pahtway IDs and elements are Entrez gene IDs.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ \code{\link{getCategoryList}}}
\examples{
\dontrun{a <- getcaBIOPATHList('1647')}
\dontrun{length(a)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
