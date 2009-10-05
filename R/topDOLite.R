`topDOLite` <-
function(x, catTerm=TRUE, ...) {
	if (length(grep('DOLite', x@categoryType)) < 1) stop(paste('The input geneAnswers categoryType is not DOLite but ', x@categoryType, '. stop function!'))
	if (catTerm) rownames(x@enrichmentInfo) <- getDOLiteTerms(rownames(x@enrichmentInfo))
	return(topCategory(x, ...))	
}

