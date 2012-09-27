`topPATH` <-
function(x, catTerm=TRUE, keepID=TRUE, ...) {
	if (length(grep('KEGG', x@categoryType)) < 1) stop(paste('The input geneAnswers categoryType is not KEGG but ', x@categoryType, '. stop function!'))
	if (catTerm) {
		if (keepID) rownames(x@enrichmentInfo) <- paste(unlist(getPATHTerms(rownames(x@enrichmentInfo))), '::', rownames(x@enrichmentInfo), sep='')
		else rownames(x@enrichmentInfo) <- unlist(getPATHTerms(rownames(x@enrichmentInfo)))
	}
	return(topCategory(x, ...))
}

