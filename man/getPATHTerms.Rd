\name{getPATHTerms}
\alias{getPATHTerms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get Pathway names of given KEGG IDs }
\description{
  Function to convert given KEGG IDs to Pathway names.
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
\references{ ~put references to the literature/web site here ~ }
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
getPATHTerms(c('04916', '05221'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
