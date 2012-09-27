`topPATHGenes` <-
function(x, catTerm=TRUE, keepID=TRUE, geneSymbol=TRUE, ...) {
	if (length(grep('KEGG', x@categoryType)) < 1) stop(paste('The input geneAnswers categoryType is not KEGG but ', x@categoryType, '. stop function!'))
	if (catTerm) {
		if (keepID) {
			rownames(x@enrichmentInfo) <- paste(unlist(getPATHTerms(rownames(x@enrichmentInfo))), '::', rownames(x@enrichmentInfo), sep='') 
			names(x@genesInCategory) <- paste(unlist(getPATHTerms(names(x@genesInCategory))), '::', names(x@genesInCategory), sep='') 
		} else {
			rownames(x@enrichmentInfo) <- unlist(getPATHTerms(rownames(x@enrichmentInfo)))
			names(x@genesInCategory) <- unlist(getPATHTerms(names(x@genesInCategory)))
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

