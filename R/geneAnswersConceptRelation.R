`geneAnswersConceptRelation` <- 
function(x, showCats=c(1:5), conceptsIDs=NULL, directed=TRUE, direction=c('down', 'up', 'both'), catTerm=TRUE, catID=FALSE, nameLength='all', ...) {
	direction <- match.arg(direction)
	enrichM <- as.matrix(getEnrichmentInfo(x))
	if (is.numeric(showCats)) {
		if (!(all(showCats %in% c(1:dim(enrichM)[1])))) print('Some specified categories might not be statistical significant! Only show significant categories.')
		showCats <- intersect(showCats, c(1:dim(enrichM)[1]))
		graphIDs <- rownames(enrichM)[showCats] 
	} else {
		if (is.character(showCats)) {
			showCats <- intersect(showCats, rownames(enrichM))
			if (length(showCats) < 1) stop('specified categories can not be recognized!') 
			graphIDs <- showCats
		} else stop('specified categories can not be recognized!')
	}
	if (is.null(conceptsIDs)) {
		conceptsIDs <- cbind(rownames(enrichM), -log2(enrichM[,dim(enrichM)[2]]), enrichM[,'genes in Category'])
	} else {
		conceptsIDs <- as.matrix(conceptsIDs)
	}
	return(getConnectedGraph(graphIDs, filterGraphIDs=conceptsIDs, idType=getCategoryType(x), directed=directed, direction=direction, readable=catTerm, nameLength=nameLength, addID=catID, ...))
}