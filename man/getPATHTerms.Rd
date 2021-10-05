\name{getPATHTerms}
\alias{getPATHTerms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get Pathway names of given KEGG IDs }
\description{
  Function to map given KEGG IDs to Pathway names.
}
\usage{
getPATHTerms(pathIDs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pathIDs}{ a KEGG IDs vector }
}
\value{
  return a KEGG pathway terms of given KEGG IDs.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{\code{\link{getCategoryTerms}}}
\examples{
getPATHTerms(c('hsa04916', 'hsa05221'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
