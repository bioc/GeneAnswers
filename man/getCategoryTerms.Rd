\name{getCategoryTerms}
\alias{getCategoryTerms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mapping Category IDs to Terms }
\description{
  Function to map category IDs to category terms.
}
\usage{
getCategoryTerms(catIDs, catType, strict = FALSE, missing=c('name', 'keep', 'remove'), nameLength='all', addID=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{catIDs}{ a character vector containing category IDs }
  \item{catType}{ type of category }
  \item{strict}{ logic value to stop conversion if NA is introduced. }
  \item{missing}{ type of handling NA mapping. }
  \item{nameLength}{ show how many first letters for long term names, 'all' for full name }
  \item{addID}{ logic, add term IDs following term names or not}
}
\details{
  The current version only supports 'GO', 'DOLITE', 'KEGG', 'REACTOME.PATH' and 'CABIO.PATH'.
  There are three types of parameters for variable 'missing'. 'name' means the NA mapping values are replaced by their names. 'keep' means all of NA values are kept. 'remove' means all of NA values are removed.
}
\value{
  return category terms of given category IDs.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\examples{
getCategoryTerms(c("04640", "05221", "05215"), catType='KEGG')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
