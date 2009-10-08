\name{mouseGeneInput}
\alias{mouseGeneInput}
\docType{data}
\title{ Example mouse gene data}
\description{
  An example of a group of mouse gene data.
}
\usage{data(mouseGeneInput)}
\format{
  A data frame with 71 observations. Column names are "Symbol", "foldChange" and "pValue". Row names are gene Entrez IDs. For general usage, row names of geneInput could be anything. 
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
data(mouseGeneInput)
mouseGeneInput[1:10,] 
}
\keyword{datasets}
