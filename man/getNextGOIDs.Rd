\name{getNextGOIDs}
\alias{getNextGOIDs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ retrieve parents or children GO IDs for given GO IDs }
\description{
  A function to retrieve parents or children GO IDs for given IDs with specified filtered IDs.
}
\usage{
getNextGOIDs(GOIDs, GOType=c('GO', 'GO.BP', 'GO.CC', 'GO.MF'), remove=TRUE, filterGOIDs=NULL, UP=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GOIDs}{ a character GO ID vector }
  \item{GOType}{ type of GO IDs, 'GO', 'GO.BP', 'GO.CC' and 'GO.MF'  }
  \item{remove}{ logic, remove the empty GOIDs in the return values}
  \item{filterGOIDs}{ a chacater vector for filtered GO IDs}
  \item{UP}{ logic, determine search Parents or Children.}
}
\details{
 filterGraphIDs is used to only keep nodes in filterGraphIDs.  
}
\value{
  return a GO IDs list representing a network.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\examples{
getNextGOIDs(c('GO:0050794','GO:0034960'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
