`topREACTOME.PATHGenes` <-
function(x, catTerm=TRUE, keepID=TRUE, geneSymbol=TRUE, ...) {
	if (length(grep('REACTOME.PATH', x@categoryType)) < 1) stop(paste('The input geneAnswers categoryType is not REACTOME.PATH but ', x@categoryType, '. stop function!'))
	if (catTerm) {
		if (keepID) {
			rownames(x@enrichmentInfo) <- paste(unlist(getREACTOMEPATHTerms(rownames(x@enrichmentInfo))), '::', rownames(x@enrichmentInfo), sep='') 
			names(x@genesInCategory) <- paste(unlist(getREACTOMEPATHTerms(names(x@genesInCategory))), '::', names(x@genesInCategory), sep='')
		} else {
			rownames(x@enrichmentInfo) <- getREACTOMEPATHTerms(rownames(x@enrichmentInfo))
			names(x@genesInCategory) <- getREACTOMEPATHTerms(names(x@genesInCategory))
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

