\name{getREACTOMEPATHList}
\alias{getREACTOMEPATHList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Retrieve REACTOME path categories containing given genes }
\description{
  Function to retrieve REACTOME path_db IDs containing given genes.
}
\usage{
getREACTOMEPATHList(geneVector, lib)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneVector}{ an Entrez gene IDs vector }
  \item{lib}{ annotation library to be used to retrieve REACTOME path_db IDs IDs. }
}
\details{
  The current version only supports Bioconductor team maintained annotation libraries, like 'org.Bt.eg.db', 'org.Ce.eg.db', 'org.Cf.eg.edu', 'org.Dm.eg.db', 'org.Dr.eg.db', 'org.EcK12.eg.db', 'org.EcSakai.eg.db', 'org.Gg.eg.db', 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db' and 'org.Ss.eg.db'.
  If the REACTOME service is not available, the function will stop.
}
\value{
  return a REACTOME genes ID list, names of the list are REACTOME path IDs IDs and elements are gene IDs.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ \code{\link{getCategoryList}}}
\examples{
\dontrun{a <- getREACTOMEPATHList(c('8772', '1017'), 'org.Hs.eg.db')}
\dontrun{length(a)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
