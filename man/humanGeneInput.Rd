\name{humanGeneInput}
\alias{humanGeneInput}
\docType{data}
\title{ Example human gene data}
\description{
  An example of a group of human gene data.
}
\usage{data(humanGeneInput)}
\format{
  A data frame with 86 observations. Column names are "Symbol", "foldChange" and "pValue". Row names are gene Entrez IDs.
}
\details{
  Fold change could be negative, which means the ratio of treatment to control is less than 1 and the value is reciprocal of general fold change.
}
\source{
  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
  ~~ possibly secondary sources and usages ~~
}
\examples{
data(humanGeneInput)
humanGeneInput[1:10,]
}
\keyword{datasets}
