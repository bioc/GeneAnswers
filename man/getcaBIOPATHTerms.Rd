\name{getcaBIOPATHTerms}
\alias{getcaBIOPATHTerms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get Pathway names of given REACTOME PATH_DB IDs }
\description{
  Function to map given caBIO pathway IDs to Pathway names. This function is not supported starting this version due to the termination of the caBig project.
}
\usage{
getcaBIOPATHTerms(caBIOPATHIDs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{caBIOPATHIDs}{ a caBIO pathway IDs vector }
}
\details{
  caBIO(Cancer Bioinformatics Infrastructure Objects, https://cabig.nci.nih.gov/tools/cabio) integrates three pathway databases from NCI-Nature curated, Biocarta and Reactome. Therefore, terms
could be same from different databases and the source library is added the end of each term. 
}
\value{
  return the caBIO pathway terms of given caBIO pathway IDs.  
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\examples{
\dontrun{getcaBIOPATHTerms(c('7622', '289', '7173'))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
