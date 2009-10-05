`geneAnswersConcepts` <- 
function(x, centroidSize=c('geneNum', 'pvalue', 'foldChange', 'oddsRatio', 'correctedPvalue'), output=c('fixed','interactive'), showCats=c(1:5), catTerm=FALSE, geneSymbol=FALSE) {
	centroidSize <- match.arg(centroidSize)
	if ((centroidSize == 'correctedPvalue') & !('fdr p value' %in% colnames(x@enrichmentInfo))) stop('input geneAnswer class does not contain fdr p value!!!')
	inputX <- geneAnswersReadable(x, catTerm=catTerm, geneSymbol=geneSymbol)
	orderby <- switch(centroidSize,
		'geneNum'= c('genes in Category', 'TRUE', 'Normal'),
		'pvalue' = c('p value', 'FALSE', '-Log10'),
		'foldChange' = c('fold of overrepresents', 'TRUE', 'Normal'),
		'oddsRatio' = c('odds ratio', 'FALSE', '-Log'), 
		'correctedPvalue' = c('fdr p value', 'FALSE', '-Log10'))
	inputX@enrichmentInfo <- inputX@enrichmentInfo[order(inputX@enrichmentInfo[, orderby[1]], decreasing=as.logical(orderby[2])), ]
	centroidSize <- inputX@enrichmentInfo[, orderby[1]]	
	names(centroidSize) <- rownames(inputX@enrichmentInfo)
	if (is.numeric(showCats)) showCats <- intersect(showCats, c(1:dim(x@enrichmentInfo)[1]))
    else {
		if (is.character(showCats)) {
			showCats <- intersect(showCats, rownames(x@enrichmentInfo))
			if (length(showCats) < 1) stop('specified categories can not be recognized!')
		}
		else stop('specified categories can not be recognized!')
	}
	temp <- centroidSize[showCats]
	scaledTemp <- switch(orderby[3],
		'Normal'=temp,
		'-Log'=-log(temp),
		'-Log10'= -log10(temp))
	categoryNet(inputX@genesInCategory[showCats], centroidSize=scaledTemp, output=output)
	return(invisible(inputX)) 
}