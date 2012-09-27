\name{getTotalGeneNumber}
\alias{getTotalGeneNumber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Obtain the total number of genes in the given annotation library }
\description{
  A function to Obtain the total number of genes in the given annotation library.
}
\usage{
getTotalGeneNumber(categoryType=c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'DOLITE', 'KEGG', 'REACTOME.PATH', 'CABIO.PATH'), known=TRUE, annotationLib=c('org.Ag.eg.db', 'org.Bt.eg.db', 'org.Ce.eg.db', 'org.Cf.eg.db', 'org.Dm.eg.db', 
							'org.Dr.eg.db', 'org.EcK12.eg.db', 'org.EcSakai.eg.db', 'org.Gg.eg.db', 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Mmu.eg.db', 'org.Pt.eg.db', 'org.Rn.eg.db', 
							'org.Ss.eg.db', 'org.Xl.eg.db', 'org.At.tair.db', 'org.Pf.plasmo.db', 'org.Sc.sgd.db'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{categoryType}{ name of given annotation category or NULL for user provided annotation list. }
  \item{known}{logic, specify only known annotation gene enrichment test. }
  \item{annotationLib}{ name of given annotation library file or user provided annotation list. }
}
\details{
  categoryType could be one of "GO", "GO.BP", "GO.CC", "GO.MF", "DOLITE", "KEGG", "REACTOME.PATH" and "CABIO.PATH".

  annotationLib could be one of "org.Ag.eg.db", "org.Bt.eg.db", "org.Ce.eg.db", "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db", "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", "org.Mmu.eg.db", "org.Pt.eg.db", "org.Rn.eg.db", "org.Ss.eg.db", "org.Xl.eg.db", "org.At.tair.db", "org.Pf.plasmo.db" and "org.Sc.sgd.db". However, if caegoryType is set to "REACTOME.PATH", only 'org.At.tair.db'(516), 'org.Ce.eg.db'(627), 'org.Dm.eg.db'(686), 'org.EcK12.eg.db'(185), 'org.EcSakai.eg.db'(185), 'org.Gg.eg.db'(840), 'org.Hs.eg.db'(1019), 'org.Mm.eg.db'(900), 'org.Pf.plasmo.db'(308), 'org.Rn.eg.db'(883) and 'org.Sc.sgd.db'(473) are available. Since DOLITE is designed for human being, currently only 4051 genes are annotated in Disease Ontology. Other species could be mapped to homologous genes by \code{\link{getHomoGeneIDs}}.

  If known is set to TRUE, the enrichment test only considers the genes with annotation. If FALSE, the total number of genes in that species will be returned.
}
\value{
  A number of total genes.
}
\references{ Feng, G., Du, P., Krett, N., Tessel, M., Rosen, S., Kibbe, W.A. and Lin, S.M., 'A collection of bioconductor methods to visualize gene-list annotations', BMC Research Notes 2010, 3:10}
\author{ Gang Feng,  Pan Du and Simon Lin }
\seealso{ \code{\link{geneAnswersBuilder}}} 
\examples{
getTotalGeneNumber(categoryType='GO.CC', annotationLib='org.Hs.eg.db')}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
