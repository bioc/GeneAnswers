`geneAnswersHeatmap` <-
function (x, showCats=c(1:5), catTerm=FALSE, geneSymbol=FALSE, ...) {
	if (is.null(x@genesInCategory[showCats])) stop('specified categories can not be found in x@genesInCategory!')
	if (is.null(x@geneExprProfile)) stop('Gene expression file is NULL!')
	if (is.numeric(showCats)) newList <- x@genesInCategory[intersect(showCats, c(1:length(x@genesInCategory)))]
    else {
		if (is.character(showCats)) newList <- x@genesInCategory[names(x@genesInCategory) %in% showCats]
		else stop('specified categories can not be recognized!')
	}
	newDataMatrix <- as.matrix(x@geneExprProfile) 
	if (geneSymbol) {
		newList <- lapply(x@genesInCategory[showCats], getSymbols, x@annLib)
		rownames(newDataMatrix) <- getSymbols(rownames(x@geneExprProfile), x@annLib)
	}
	if (is.numeric(showCats)) showCats <- intersect(showCats, c(1:dim(x@enrichmentInfo)[1]))
	else {
		if (is.character(showCats)) showCats <- intersect(showCats, rownames(x@enrichmentInfo))
		else stop('specified categories can not be recognized!')
	}
	if (catTerm) names(newList) <- getCategoryTerms(names(x@genesInCategory[showCats]), x@categoryType)
	geneAnnotationHeatmap(newList, dataMatrix = newDataMatrix, ...) 		  
}

