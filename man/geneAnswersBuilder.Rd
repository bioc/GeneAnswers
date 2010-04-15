\name{geneAnswersBuilder}
\alias{geneAnswersBuilder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Build an object of a GeneAnswers class }
\description{
  A function to build an object of a GeneAnswers class based on given information.
}
\usage{
geneAnswersBuilder(geneInput, annotationLib, categoryType = NULL, testType = c("hyperG", "none"), totalGeneNumber=NULL, geneExpressionProfile = NULL, categorySubsetIDs = NULL, pvalueT = 0.01, FDR.correction = FALSE, verbose=TRUE, sortBy=c('pvalue', 'geneNum', 'foldChange', 'oddsRatio', 'correctedPvalue', 'none'), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneInput}{ a dataframe containing gene IDs and possible values associated with given gene IDs. }
  \item{annotationLib}{ name of given annotation library file or user provided annotation list. }
  \item{categoryType}{ name of given annotation category or NULL for user provided annotation list. }
  \item{testType}{ name of enrichment test. }
  \item{totalGeneNumber}{ number of total genes to perform hypergeometric test. }
  \item{geneExpressionProfile}{ data frame containing gene expression file or NULL. }
  \item{categorySubsetIDs}{ a character vector of user-specified subset of categories to be tested. }
  \item{pvalueT}{ p-value threshold of the enrichment test. }
  \item{FDR.correction}{ logical indicating if FDR correction of the enrichment test p-value is performed or not. }
  \item{verbose}{ logical, display current building stage. }
  \item{sortBy}{ sorted type } 
  \item{\dots}{ additional arguments passed on to the corresponding enrichment test. }
}
\details{
  As the input of geneAnswersBuilder, geneInput could be a character vector (Gene Entrez ID vector), a matrix or a dataframe. For the matrix and dataframe, the first column is for Gene Entrez IDs, while other columns could be any interested values that could be used to represent gene expression direction for generating concepts-genes network. Rownames are not necessary. 
  
  annotationLib could be Disease Ontology library, Entrez annotation libraries for a specie, such as 'org.Hs.eg.db'. Current version supports 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db' and 'org.Dm.eg.db'. User can also use own annotation library. User's annotation library should be a list. Each element in this list is a vector of genes for a user-specified category. Names of this annotation list are categories' names.

  categoryType could be "GO", "GO.BP", "GO.CC", "GO.MF", "DOLite", "KEGG". "GO.BP" only test biological process Gene Ontology terms, "GO.CC" for cellular components, "GO.MF" for molecular functions, and "GO" for all of these three categories. For user provided annotation library, it should be NULL in most cases.
  
  totalGeneNumber could be NULL if annotationLib is one of 'org.Hs.eg.db'(45384), 'org.Mm.eg.db'(61498), 'org.Rn.eg.db'(37536) and 'org.Dm.eg.db'(22606). If user has own annotationLib, totalGeneNumber should be an integer, or one of 'human', 'mouse', 'rat' and 'fly'. geneAnswersBuilder will automatically assign the corresponding value to totalGeneNumber.  

  sortBy could be one of "geneNum", "pvalue", "foldChange", "oddsRatio", "correctedPvalue" and "none". Default value is 'pvalue'.
}
\value{
  A GeneAnswers class containing geneInput, entrichmentInfo, etc.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng,  Pan Du and Simon Lin }
\seealso{ \code{\link{geneAnswersSort}}, \code{\link{geneAnswersConceptNet}}} 
\examples{
data('humanExpr')
data('humanGeneInput') 
x <- geneAnswersBuilder(humanGeneInput, 'org.Hs.eg.db', categoryType='GO.BP', testType='h', pvalueT=0.1, FDR.correct=TRUE, geneExpressionProfile=humanExpr)
class(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
