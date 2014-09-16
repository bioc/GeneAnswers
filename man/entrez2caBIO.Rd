\name{entrez2caBIO}
\alias{entrez2caBIO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ map Entrez gene IDs to caBIO gene IDs}
\description{
  Function to map the given Entrez gene IDs to the caBIO gene IDs. This function is not supported starting this version due to the termination of the caBig project.
}
\usage{
entrez2caBIO(lls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lls}{ an Entrez gene IDs vector }
}
\value{
  return a caBIO genes ID list, names of the list are Entrez gene IDs and elements are caBIO gene IDs.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\examples{
\dontrun{entrez2caBIO(c('1647', '596'))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
