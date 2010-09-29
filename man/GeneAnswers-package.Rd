\name{GeneAnswers-package}
\alias{GeneAnswers-package}
\alias{GeneAnswers}
\docType{package}
\title{
Integrated Interpretation of Genes
}
\description{
GeneAnswers provide an integrated tool for biological or medical interpretation of the given one or more groups of genes by means of statistical test. 
}
\details{
\tabular{ll}{
Package: \tab GeneAnswers\cr
Type: \tab Package\cr
Version: \tab 1.6.0\cr
Date: \tab 2010-10-14\cr
License: \tab LGPL version 2 or newer \cr
}
}
\author{
Gang Feng, Pan Du, Tian Xia, Warren Kibbe and Simon Lin


Maintainer: Gang Feng <g-feng@northwestern.edu> and Pan Du <dupan@northwestern.edu>
}
\references{
1.	Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10

2. Du, P., Feng, G., Flatow, J., Song, J., Holko, M., Kibbe, W.A. and Lin, S.M., (2009) 'From disease ontology to disease-ontology lite: statistical methods to adapt a general-purpose ontology for the test of gene-ontology associations', Bioinformatics 25(12):i63-8

3. Osborne, J.D., Flatow, J., Holko, M., Lin, S.M., Kibbe, W.A., Zhu, L.J., Danila, M.I., Feng, G. and Chisholm, R.L., Annotating the human genome with Disease Ontology. BMC Genomics. 2009 Jul 7;10 Suppl 1:S6.
}
\examples{
data('humanExpr')
data('humanGeneInput')
x <- geneAnswersBuilder(humanGeneInput, 'org.Hs.eg.db', categoryType='GO.BP', testType='hyperG', pvalueT=0.01, FDR.correct=TRUE, geneExpressionProfile=humanExpr)
class(x)
}
\keyword{ package }