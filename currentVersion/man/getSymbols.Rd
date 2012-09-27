\name{getSymbols}
\alias{getSymbols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert entrez gene IDs to gene symbols }
\description{
  function to convert given entrez gene IDs to gene symbols. 
}
\usage{
getSymbols(geneIDs, data, strict = FALSE, missing=c('name', 'keep', 'remove'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneIDs}{ an Entrez gene IDs vector }
  \item{data}{ annotation library }
  \item{strict}{ logic value to stop conversion if NA is introduced. }
  \item{missing}{ type of handling NA mapping. }
}

\value{
  return a gene symbols vector of given gene IDs.
  There are three types of parameters for variable 'missing'. 'name' means the NA mapping values are replaced by their names. 'keep' means all of NA values are kept. 'remove' means all of NA values are removed.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\examples{
require('org.Mm.eg.db')
getSymbols(c('11651', '11836'), 'org.Mm.eg.db')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
