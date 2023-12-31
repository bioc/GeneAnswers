`topDOLITE` <-
function(x, catTerm=TRUE, keepID=TRUE, ...) {
	if (length(grep('DOLITE', x@categoryType)) < 1) stop(paste('The input geneAnswers categoryType is not DOLite but ', x@categoryType, '. stop function!'))
	if (catTerm) {
		if (keepID) rownames(x@enrichmentInfo) <- paste(getDOLiteTerms(rownames(x@enrichmentInfo)), '::', rownames(x@enrichmentInfo), sep='')
		else rownames(x@enrichmentInfo) <- getDOLiteTerms(rownames(x@enrichmentInfo))
	}
	return(topCategory(x, ...))	
}

