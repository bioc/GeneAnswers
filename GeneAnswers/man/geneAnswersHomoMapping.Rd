\name{geneAnswersHomoMapping}
\alias{geneAnswersHomoMapping}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mapping homogenes for a GeneAnswers instance }
\description{
  A function to mapping homogenes in all of slots of a GeneAnswer instance
}
\usage{
geneAnswersHomoMapping(x, species = c("human", "rat", "mouse", "fly"), speciesL = c("human", "rat", "mouse", "fly"), mappingMethod = c("direct", "biomaRt", "none"), filterGenes = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a GeneAnswers instance }
  \item{species}{ species of the current genes }
  \item{speciesL}{ species of the mapped genes }
  \item{mappingMethod}{ mapping method, see details }
  \item{filterGenes}{ a gene symbol vector to filter genes}
  \item{verbose}{ logical, show current stage or not}
}
\details{
  There are two mapping methods supported by current version. "direct" only works between human and mouse because most of human gene symbols are capitalized and only the first letter is uppercase for those homogenes in mouse. Another way is by means of package "biomaRt" ,which contains more information while the network connection is necessary to access biomaRt online server.
  Since two methods are based on different mechanisms, it is highly recommended to employ same method during mapping.
  Each method might introduce more homogenes, so users can remove ones that do not belong to original genes by optional "filterGeneList". 
}
\value{
  return a mapped GeneAnswers instance
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ \code{\link{getHomoGeneIDs}}}
\examples{
example(GeneAnswers)
\dontrun{geneAnswersHomoMapping(x, species='human', speciesL='mouse', mappingMethod='direct')}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
