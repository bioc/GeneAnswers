\name{getGOList}
\alias{getGOList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get GO list of given genes }
\description{
  Retrieve GO IDs based on given gene IDs.
}
\usage{
getGOList(geneVector, lib, GOCat = c("ALL", "BP", "CC", "MF"), level = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneVector}{ a character vector containing entrez IDs }
  \item{lib}{ annotation library }
  \item{GOCat}{ type of Gene Ontology }
  \item{level}{ positive integer to specify how many levels GO IDs will be removed. }
}
\details{
  User can specify which subtype of GO can be kept. "ALL" means all of subtypes are kept.
  Gene Ontology is a tree-like structure. Level can be used to remove top noncritical GO IDs.
}
\value{
  return a GO list, whose names are GO IDs. Elements are gene entrez IDs belonging to the corresponding GO categories.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
a <- getGOList(c('56458', '16590'), 'org.Mm.eg.db', GOCat='BP', level=2)
length(a)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
