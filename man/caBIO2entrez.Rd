\name{caBIO2entrez}
\alias{caBIO2entrez}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ map caBIO gene IDs to Entrez gene IDs}
\description{
  Function to map the given caBIO gene IDs to the Entrez gene IDs. This function is not supported starting this version due to the termination of the caBig project.
}
\usage{
caBIO2entrez(caBIOIds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{caBIOIds}{ an caBIOIds gene IDs vector }
}
\value{
  return a Entrez genes ID list, names of the list are the given caBIO gene IDs and elements are Entrez gene IDs.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\examples{
\dontrun{caBIO2entrez(c('2933', '7326'))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
