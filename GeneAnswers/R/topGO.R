`topGO` <-
function(x, catTerm=TRUE, keepID=TRUE, ...) {
	if (length(grep('GO', x@categoryType)) < 1) stop(paste('The input geneAnswers categoryType is not GO but ', x@categoryType, '. stop function!'))
	if (catTerm) {
		if (keepID) rownames(x@enrichmentInfo) <- paste(unlist(.getGOTerms(rownames(x@enrichmentInfo))), '::', rownames(x@enrichmentInfo), sep='')
		else rownames(x@enrichmentInfo) <- unlist(.getGOTerms(rownames(x@enrichmentInfo)))
	} 
	return(topCategory(x, ...))
}

