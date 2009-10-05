`topGO` <-
function(x, catTerm=TRUE, ...) {
	if (length(grep('GO', x@categoryType)) < 1) stop(paste('The input geneAnswers categoryType is not GO but ', x@categoryType, '. stop function!'))
	if (catTerm) rownames(x@enrichmentInfo) <- unlist(getGOTerm(rownames(x@enrichmentInfo)))
	return(topCategory(x, ...))
}

