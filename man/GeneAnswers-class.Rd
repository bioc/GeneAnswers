\name{GeneAnswers-class}
\docType{class}
\alias{class:GeneAnswers}
\alias{GeneAnswers-class}

\alias{getGeneInput,GeneAnswers-method}
\alias{getGeneInput}
\alias{getTestType,GeneAnswers-method}
\alias{getTestType}
\alias{getPValueT,GeneAnswers-method}
\alias{getPValueT}
\alias{getGenesInCategory,GeneAnswers-method}
\alias{getGenesInCategory}
\alias{getGeneExprProfile,GeneAnswers-method}
\alias{getGeneExprProfile}
\alias{getAnnLib,GeneAnswers-method}
\alias{getAnnLib}
\alias{getCategoryType,GeneAnswers-method}
\alias{getCategoryType}
\alias{getEnrichmentInfo,GeneAnswers-method}
\alias{getEnrichmentInfo}
\alias{setGeneInput,GeneAnswers-method}
\alias{setGeneInput}
\alias{setTestType,GeneAnswers-method}
\alias{setTestType}
\alias{setPValueT,GeneAnswers-method}
\alias{setPValueT}
\alias{setGenesInCategory,GeneAnswers-method}
\alias{setGenesInCategory}
\alias{setGeneExprProfile,GeneAnswers-method}
\alias{setGeneExprProfile}
\alias{setAnnLib,GeneAnswers-method}
\alias{setAnnLib}
\alias{setCategoryType,GeneAnswers-method}
\alias{setCategoryType}
\alias{setEnrichmentInfo,GeneAnswers-method}
\alias{setEnrichmentInfo}

\alias{summary,GeneAnswers-method}
\alias{summary}
\alias{show,GeneAnswers-method}
\alias{show}




\title{Class GeneAnswers: contain and describe the relationship between given gene data and specified category}
\description{This is a class representation of the relationship between given gene data and specified category.}

\section{Creating Objects}{
  Objects can be created using the function \code{\link{geneAnswersBuilder}}.
}
\section{Slots}{
   Slot specific to \code{GeneAnswers}:
   \describe{
   		\item{\code{geneInput}:}{ a data frame containing gene Entrez IDs with or without any values. Current version only supports gene Entrez IDs. The values could be foldChange, p value, or other values. These data can be used for concept-gene network. Genes with positive values will be represented as red nodes, while negative value genes are green nodes.}
   		\item{\code{testType}:}{ statistical test method. Current version supports hypergeometric test to test relationship between genes and specified categories.}
   		\item{\code{pvalueT}:}{ the cutoff value of statistical test. Any categories will not be reported if the p value is more than the cutoff.}
         \item{\code{genesInCategory}:}{ a list containing genes belonging to categories. The names of the list are categories.}
         \item{\code{geneExprProfile}:}{ a data frame to store gene expression data. If not available, it could be NULL.}
         \item{\code{annLib}:}{ annotation database used for statistical test.}
         \item{\code{categoryType}:}{functional or medical category used for statistical test.}
         \item{\code{enrichmentInfo}:}{a data frame containing filtered categories with statistical results by specified pvalueT.}
   }
}

\section{Methods}{
  \bold{Class-specific methods:}
  \describe{
     \item{\code{getGeneInput(GeneAnswers)}:}{Access the geneInput slot of \code{GeneAnswers} object.}
     \item{\code{getTestType(GeneAnswers)}:}{Access the testType slot of \code{GeneAnswers} object.}
     \item{\code{getPValueT(GeneAnswers)}:}{Access the pvaluteT slot of \code{GeneAnswers} object.}
     \item{\code{getGenesInCategory(GeneAnswers)}:}{Access the genesInCategory slot of \code{GeneAnswers} object.}
     \item{\code{getGeneExprProfile(GeneAnswers)}:}{Access the geneExprProfile slot of \code{GeneAnswers} object.}
     \item{\code{getAnnLib(GeneAnswers)}:}{Access the annLib slot of \code{GeneAnswers} object.}
     \item{\code{getCategoryType(GeneAnswers)}:}{Access the categoryType slot of \code{GeneAnswers} object.}
     \item{\code{getEnrichmentInfo(GeneAnswers)}:}{Access the enrichmentInfo slot of \code{GeneAnswers} object.}
     \item{\code{setGeneInput(GeneAnswers, geneInput)}:}{Assign the geneInput slot of \code{GeneAnswers} object.}
     \item{\code{setTestType(GeneAnswers, type=c('hyperG', 'none'))}:}{Assign the testType slot of \code{GeneAnswers} object.}
     \item{\code{setPValueT(GeneAnswers, pvalueT)}:}{Assign the pvaluteT slot of \code{GeneAnswers} object.}
     \item{\code{setGenesInCategory(GeneAnswers, genesInCategory)}:}{Assign the genesInCategory slot of \code{GeneAnswers} object.}
     \item{\code{setGeneExprProfile(GeneAnswers, geneExprProfile)}:}{Assign the geneExprProfile slot of \code{GeneAnswers} object.}
     \item{\code{setAnnLib(GeneAnswers, annLib)}:}{Assign the annLib slot of \code{GeneAnswers} object.}
     \item{\code{setCategoryType(GeneAnswers, type=c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'DOLITE', 'KEGG', 'REACTOME.PATH', 'CABIO.PATH', 'User defiend'))}:}{Assign the categoryType slot of \code{GeneAnswers} object.}
     \item{\code{setEnrichmentInfo(GeneAnswers, enrichmentInfo)}:}{Assign the enrichmentInfo slot of \code{GeneAnswers} object.}
     \item{\code{summary(GeneAnswers)}:}{Briefly summarize the information of \code{GeneAnswers} object and show contents of \code{GeneAnswers} object.}
     \item{\code{show(GeneAnswers)}:}{Briefly show contents of \code{GeneAnswers} object.}
  }
}

\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng, Pan Du and Simon Lin }
\seealso{ 
  \code{\link{geneAnswersBuilder}}
}
\examples{
data('humanExpr')
data('humanGeneInput')
x <- geneAnswersBuilder(humanGeneInput, 'org.Hs.eg.db', categoryType='GO.BP', testType='hyperG', pvalueT=0.01, FDR.correct=TRUE, geneExpressionProfile=humanExpr)
class(x)
}

\keyword{classes}
