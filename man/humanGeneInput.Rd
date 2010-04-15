\name{humanGeneInput}
\alias{humanGeneInput}
\docType{data}
\title{ Example human gene data}
\description{
  An example of a group of human gene data.
}
\usage{data(humanGeneInput)}
\format{
  A data frame with 86 observations. Column names are "Symbol", "foldChange" and "pValue". Row names are gene Entrez IDs. For general usage, row names of geneInput could be anything.
}
\details{
  Fold change could be negative, which means the ratio of treatment to control is less than 1 and the value is reciprocal of general fold change.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}\examples{
data(humanGeneInput)
humanGeneInput[1:10,]
}
\keyword{datasets}
