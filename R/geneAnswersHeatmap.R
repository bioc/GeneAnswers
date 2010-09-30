`geneAnswersHeatmap` <-
function (x, showCats=c(1:5), catTerm=FALSE, geneSymbol=FALSE, catID=FALSE, nameLength='all', showAllGenes=FALSE, ...) {
	if (is.null(x@genesInCategory[showCats])) stop('specified categories can not be found in x@genesInCategory!')
	#if (is.null(x@geneExprProfile)) stop('Gene expression file is NULL!')
	if (is.numeric(showCats)) {
		if (!(all(showCats %in% c(1:dim(x@enrichmentInfo)[1])))) print('Some specified categories might not be statistical significant! Only show significant categories.')
		showCats <- intersect(showCats, c(1:dim(x@enrichmentInfo)[1]))
	}else {
		if (is.character(showCats)) {
			showCats <- intersect(showCats, rownames(x@enrichmentInfo))
			if (length(showCats) < 1) stop('specified categories can not be recognized!')
		} else stop('specified categories can not be recognized!')
	}
	newList <- x@genesInCategory[showCats]
	
	
	if (is.null(x@geneExprProfile)) newDataMatrix <- NULL
	else {
		newDataMatrix <- matrix(as.numeric(as.matrix(x@geneExprProfile[,2:(dim(x@geneExprProfile)[2])])), ncol=(dim(x@geneExprProfile)[2]-1))
		colnames(newDataMatrix) <- colnames(x@geneExprProfile)[2:dim(x@geneExprProfile)[2]]
		rownames(newDataMatrix) <- x@geneExprProfile[,1]
	}
	if (catTerm) {
		if (x@categoryType %in% c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'DOLITE', 'KEGG', 'REACTOME.PATH', 'CABIO.PATH')) {
	  		names(newList) <- getCategoryTerms(names(newList), x@categoryType, missing='name', nameLength=nameLength, addID=catID)
		}else {
		 	print('Slot categoryType is not recognized! No mapping ...')
		}
	}
	if (showAllGenes) newList <- c(newList, 'All Genes'=list(x@geneInput[,1])) 
	if (geneSymbol) {
		newList <- lapply(newList, getSymbols, x@annLib, missing='name')
		if (!is.null(newDataMatrix)) rownames(newDataMatrix) <- getSymbols(rownames(newDataMatrix), x@annLib, missing='name')
	}
	geneAnnotationHeatmap(newList, dataMatrix = newDataMatrix, ...) 		  
}

