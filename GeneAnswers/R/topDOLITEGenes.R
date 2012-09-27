`topDOLITEGenes` <-
function(x, catTerm=TRUE, keepID=TRUE, geneSymbol=TRUE, ...) {
	if (length(grep('DOLITE', x@categoryType)) < 1) stop(paste('The input geneAnswers categoryType is not DOLite but ', x@categoryType, '. stop function!'))
	if (catTerm) {
		if (keepID) {
			rownames(x@enrichmentInfo) <- paste(unlist(getDOLiteTerms(rownames(x@enrichmentInfo))), '::', rownames(x@enrichmentInfo), sep='') 
			names(x@genesInCategory) <- paste(unlist(getDOLiteTerms(names(x@genesInCategory))), '::', names(x@genesInCategory), sep='')
		} else {
			rownames(x@enrichmentInfo) <- getDOLiteTerms(rownames(x@enrichmentInfo))
			names(x@genesInCategory) <- getDOLiteTerms(names(x@genesInCategory))
		}
	} 
	if (geneSymbol) {
		if (is.null(x@annLib)) stop('annotation library is not specified! original gene IDs will be kept.')
		else {
			x@genesInCategory <- lapply(x@genesInCategory, getSymbols, x@annLib)
			x@geneInput[,1] <- getSymbols(x@geneInput[,1], x@annLib)
		}
	}
	return(topCategoryGenes(x, ...))
}

