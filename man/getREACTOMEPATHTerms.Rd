\name{getREACTOMEPATHTerms}
\alias{getREACTOMEPATHTerms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get Pathway names of given REACTOME PATH_DB IDs }
\description{
  Function to convert given REACTOME PATH_DB IDs to Pathway names.
}
\usage{
getREACTOMEPATHTerms(pathDBIDs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pathDBIDs}{ a REACTOME PATH_DB IDs vector }
}
\value{
  return a REACTOME pathway terms of given REACTOME PATH_DB IDs.
  If the REACTOME service is not available, the function will stop.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\examples{
\dontrun{getREACTOMEPATHTerms(c('174143', '453274'))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
