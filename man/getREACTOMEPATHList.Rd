\name{getREACTOMEPATHList}
\alias{getREACTOMEPATHList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Retrieve REACTOME path categories containing given genes }
\description{
  Function to retrieve REACTOME path_db IDs containing given genes.
}
\usage{
getREACTOMEPATHList(geneVector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneVector}{ an Entrez gene IDs vector }
}
\details{
  geneVector should be a vector of Entrez IDs.
}
\value{
  return a REACTOME genes ID list, names of the list are REACTOME path IDs and elements are gene IDs.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ \code{\link{getCategoryList}}}
\examples{
\dontrun{a <- getREACTOMEPATHList(c('8772', '1017'))}
\dontrun{length(a)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
