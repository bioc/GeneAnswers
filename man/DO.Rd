\name{DO}
\alias{DO}
\alias{DO.graph.gene}
\alias{DO.graph.closure.gene}
\alias{DO2gene.map}
\alias{gene2DO.map}
\alias{DO.terms}
\docType{data}
\title{
Several data objects related with DO (Disease Ontology) and its mapping to genes
}
\description{
Several data objects related with DO (Disease Ontology) and its mapping to genes
}
\usage{data(DO)}

\details{
The data file "DO.rda" includes five datasets:

	DO.graph.gene: a graphNEL object, which shows the ontology relations of DO
	
	DO.graph.closure.gene: a graphNEL object, whose edges represent the link between a DO term and its offspring ontology terms. Only the DO terms with gene mappings were included.
	
	DO2gene.map: a list show the mapping from DOIDs to genes
	
	gene2DO.map: a list show the mapping from genes to DOIDs
	
	DO.terms: a named character vector. Its names are DOIDs and elements are DO.terms

}

\examples{
data(DO)

datasets <- c("DO.graph.gene", "DO.graph.closure.gene", "DO2gene.map", "gene2DO.map", "DO.terms")
# check the existence of these datasets:
sapply(datasets, exists)

}
\keyword{datasets}
