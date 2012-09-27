\name{getCategoryList}
\alias{getCategoryList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Retrieve categories containing given genes }
\description{
  Function to retrieve specified category IDs containing given genes.
}
\usage{
getCategoryList(geneVector, lib, categoryType)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneVector}{ an Entrez gene IDs vector }
  \item{lib}{ annotation library to be used to retrieve categories terms. }
  \item{categoryType}{ type of category }
}
\details{
  The current version only supports Bioconductor team maintained annotation libraries, like 'org.Bt.eg.db', 'org.Ce.eg.db', 'org.Cf.eg.edu', 'org.Dm.eg.db', 'org.Dr.eg.db', 'org.EcK12.eg.db', 'org.EcSakai.eg.db', 'org.Gg.eg.db', 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db' and 'org.Ss.eg.db'.
}
\value{
  return a category list, names of the list are category IDs and elements are genes IDs.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\examples{
getCategoryList(c('56458', '16590'), 'org.Mm.eg.db', 'PATH')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
