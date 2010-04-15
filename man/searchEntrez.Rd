\name{searchEntrez}
\alias{searchEntrez}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Search specified information from Entrez site }
\description{
  A function to search Entrez website by one given keywords list.
}                                           
\usage{
searchEntrez(tagList, species = "human")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tagList}{ keyword list to search on Entrez. }
  \item{species}{ specie for search on Entrez. }
}
\value{
  an Entrez ID list containing all of relative genes from Entrez database.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Pan Du, Gang Feng and Simon Lin }
\examples{
tagList <- list(FSHR=c('"FSHR"', '"Follicle stimulating hormone receptor"'), apoptosis=c('apoptosis'))
\dontrun{entrezList <- searchEntrez(tagList, species='mouse')}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
