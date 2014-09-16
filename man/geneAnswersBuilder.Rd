\name{geneAnswersBuilder}
\alias{geneAnswersBuilder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Build an object of a GeneAnswers class }
\description{
  A function to build an object of a GeneAnswers class based on given information.
}
\usage{
geneAnswersBuilder(geneInput, annotationLib, categoryType = NULL, testType = c("hyperG", "none"), known=TRUE, totalGeneNumber=NULL, geneExpressionProfile = NULL, categorySubsetIDs = NULL, pvalueT = 0.01, FDR.correction = FALSE, verbose=TRUE, sortBy=c('pvalue', 'geneNum', 'foldChange', 'oddsRatio', 'correctedPvalue', 'none'), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneInput}{ a dataframe containing gene IDs and possible values associated with given gene IDs. }
  \item{annotationLib}{ name of given annotation library file or user provided annotation list. }
  \item{categoryType}{ name of given annotation category or NULL for user provided annotation list. }
  \item{testType}{ name of enrichment test. }
  \item{known}{logic, specify only known annotation gene enrichment test. }
  \item{totalGeneNumber}{ number of total genes to perform hypergeometric test. }
  \item{geneExpressionProfile}{ data frame containing gene expression file or NULL. }
  \item{categorySubsetIDs}{ a character vector of user-specified subset of categories to be tested. }
  \item{pvalueT}{ p-value threshold of the enrichment test. }
  \item{FDR.correction}{ logic, indicating if FDR correction of the enrichment test p-value is performed or not. }
  \item{verbose}{ logic, display current building stage. }
  \item{sortBy}{ sorted type } 
  \item{\dots}{ additional arguments passed to \code{\link{getGOList}}. }
}
\details{
  As the input of geneAnswersBuilder, geneInput could be a character vector (Gene Entrez ID vector), a matrix or a dataframe. For the matrix and dataframe, the first column is for Gene Entrez IDs, while other columns could be any interested values that could be used to represent gene expression direction for generating concepts-genes network. Rownames are not necessary. 
  
  annotationLib could be Disease Ontology library, Entrez annotation libraries for a specie, such as 'org.Hs.eg.db'. Current version supports "org.Ag.eg.db", "org.Bt.eg.db", "org.Ce.eg.db", "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db", "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", "org.Mmu.eg.db", "org.Pt.eg.db", "org.Rn.eg.db", "org.Ss.eg.db", "org.Xl.eg.db", "org.At.tair.db", "org.Pf.plasmo.db" and "org.Sc.sgd.db". User can also use own annotation library. User's annotation library should be a list. Each element in this list is a vector of genes for a user-specified category. Names of this annotation list are categories' names.

  categoryType could be "GO", "GO.BP", "GO.CC", "GO.MF", "DOLITE", "KEGG", "REACTOME.PATH" and "CABIO.PATH". "GO.BP" only test biological process Gene Ontology terms, "GO.CC" for cellular components, "GO.MF" for molecular functions, "GO" for all of these three categories, "KEGG" for all KEGG pathways, and "REACTOME.PATH" for all REACTOME pathways. For user provided annotation library, it should be NULL in most cases.

"caBIO.PATH", which was for NCI-Nature curated, Biocarta and REACTOME, is not supported starting this version because of the termination of caBig project.
  
  If known is set to TRUE, the enrichment test only considers the genes with annotation. If FALSE, the total number of genes in that species will be returned. If user has own annotationLib, totalGeneNumber should be an integer, or one of "anopheles", "arabidopsis", "bovine", "worm", "canine", "fly", "zebrafish", "ecolistraink12", "ecolistrainsakai", "chicken", "human", "mouse", "rhesus", "malaria", "chimp", "rat", "yeast", "pig" and "xenopus". NULL only works when "known" is set TRUE. geneAnswersBuilder will automatically assign the corresponding value to totalGeneNumber. User can get total gene numbers by \code{\link{getTotalGeneNumber}}, too. 

  sortBy could be one of "geneNum", "pvalue", "foldChange", "oddsRatio", "correctedPvalue" and "none". Default value is 'pvalue'.
}
\value{
  A GeneAnswers class containing geneInput, entrichmentInfo, etc.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng,  Pan Du and Simon Lin }
\seealso{ \code{\link{getTotalGeneNumber}}} 
\examples{
data('humanExpr')
data('humanGeneInput') 
x <- geneAnswersBuilder(humanGeneInput, 'org.Hs.eg.db', categoryType='GO.BP', testType='h', pvalueT=0.1, FDR.correct=TRUE, geneExpressionProfile=humanExpr)
class(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
