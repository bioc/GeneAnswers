`topPATH` <-
function(x, catTerm=TRUE, ...) {
	if (length(grep('KEGG', x@categoryType)) < 1) stop(paste('The input geneAnswers categoryType is not KEGG but ', x@categoryType, '. stop function!'))
	if (catTerm) rownames(x@enrichmentInfo) <- unlist(getPATHTerms(rownames(x@enrichmentInfo)))
	return(topCategory(x, ...))
}

