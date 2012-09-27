`geneAnswersConcepts` <- 
function(x, centroidSize=c('geneNum', 'pvalue', 'foldChange', 'oddsRatio', 'correctedPvalue'), output=c('fixed','interactive'), showCats=c(1:5), catTerm=FALSE, catID=FALSE) {
	centroidSize <- match.arg(centroidSize)
	if ((centroidSize == 'correctedPvalue') & !('fdr p value' %in% colnames(x@enrichmentInfo))) stop('input geneAnswer class does not contain fdr p value!!!')
	if (is.numeric(showCats)) {
		if (!(all(showCats %in% c(1:dim(x@enrichmentInfo)[1])))) print('Some specified categories might not be statistical significant! Only show significant categories.')
		showCats <- intersect(showCats, c(1:dim(x@enrichmentInfo)[1])) 
	} else {
		if (is.character(showCats)) {
			showCats <- intersect(showCats, rownames(x@enrichmentInfo))
			if (length(showCats) < 1) stop('specified categories can not be recognized!')
		} else stop('specified categories can not be recognized!')
	}
	orderby <- switch(centroidSize,
		'geneNum'= c('genes in Category', 'TRUE', 'Normal'),
		'pvalue' = c('p value', 'FALSE', '-Log10'),
		'foldChange' = c('fold of overrepresents', 'TRUE', 'Normal'),
		'oddsRatio' = c('odds ratio', 'FALSE', '-Log'), 
		'correctedPvalue' = c('fdr p value', 'FALSE', '-Log10'))
	x@enrichmentInfo <- x@enrichmentInfo[order(x@enrichmentInfo[, orderby[1]], decreasing=as.logical(orderby[2])), ]
	centroidSize <- x@enrichmentInfo[, orderby[1]]	
	names(centroidSize) <- rownames(x@enrichmentInfo)
	
	temp <- centroidSize[showCats]
	scaledTemp <- switch(orderby[3],
		'Normal'=temp,
		'-Log'=-log(temp),
		'-Log10'= -log10(temp))
    newList <- x@genesInCategory[showCats]
	if (catTerm) {
		if (x@categoryType %in% c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'DOLITE', 'KEGG', 'REACTOME.PATH', 'CABIO.PATH')) {
			if (catID) {
				names(newList) <- paste(getCategoryTerms(names(newList), x@categoryType, missing='name'), '::', names(newList), sep='')
				names(scaledTemp) <- paste(getCategoryTerms(names(scaledTemp), x@categoryType, missing='name'), '::', names(scaledTemp), sep='')
			} else {
				names(newList) <- getCategoryTerms(names(newList), x@categoryType, missing='name')
				names(scaledTemp) <- getCategoryTerms(names(scaledTemp), x@categoryType, missing='name')
			}
		} else {
			print('Slot categoryType is not recognized! No mapping ...')
		}
	}
	
	names(newList) <- paste(names(newList), '::', sapply(newList, length), sep='')
	
	categoryNet(newList, centroidSize=scaledTemp, output=output)
	return(invisible(x)) 
}