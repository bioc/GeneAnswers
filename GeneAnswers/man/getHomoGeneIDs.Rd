\name{getHomoGeneIDs}
\alias{getHomoGeneIDs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get homologous genes of given genes }
\description{
  Map given gene IDs to homologous gene IDs.
}
\usage{
getHomoGeneIDs(oriGeneIDs, species = c("human", "rat", "mouse", "yeast", "fly"), speciesL = c("human", "rat", "mouse", "yeast", "fly"), mappingMethod = c("direct", "biomaRt", "none"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oriGeneIDs}{ a given entrez gene IDs }
  \item{species}{ species of the current genes }
  \item{speciesL}{ species of the mapped genes }
  \item{mappingMethod}{ mapping method, see details }
}
\details{
  There are two mapping methods supported by current version. "direct" only works between human and mouse because most of human gene symbols are capitalized and only the first letter is uppercase for those homogenes in mouse. Another way is by means of package "biomaRt" ,which contains more information while the network connection is necessary to access biomaRt online server.
}
\value{
  return homologous gene IDs of given genes
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\examples{
getHomoGeneIDs(c('56458', '16590'), species='m', speciesL='h', mappingMethod='direct')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
