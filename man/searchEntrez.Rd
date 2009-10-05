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
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  an Entrez ID list containing all of relative genes from Entrez database.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Pan Du, Gang Feng and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
tagList <- list(FSHR=c('"FSHR"', '"Follicle stimulating hormone receptor"'), apoptosis=c('apoptosis'))
entrezList <- searchEntrez(tagList, species='mouse')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
