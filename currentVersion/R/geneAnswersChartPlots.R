`geneAnswersChartPlots` <-
function(x, chartType=c('pieChart', 'barPlot', 'all'), sortBy = c('geneNum', 'pvalue', 'foldChange', 'oddsRatio', 'correctedPvalue'), newWindow=TRUE, ...) {
	chartType <- match.arg(chartType)
	sortBy <- match.arg(sortBy)
	if ((sortBy == 'correctedPvalue') & !('fdr p value' %in% colnames(x@enrichmentInfo))) stop('input GeneAnswer instance does not contain corrected p value!!!')
	orderby <- switch(sortBy,
		'geneNum'= c('genes in Category', 'TRUE'),
		'pvalue' = c('p value', 'FALSE'),
		'foldChange' = c('fold of overrepresents', 'TRUE'),
		'oddsRatio' = c('odds ratio', 'FALSE'), 
		'correctedPvalue' = c('fdr p value', 'FALSE'))
  	y <- x
  	y@enrichmentInfo <-	x@enrichmentInfo[order(x@enrichmentInfo[, orderby[1]], decreasing=as.logical(orderby[2])), ]
  	y@genesInCategory <- x@genesInCategory[rownames(y@enrichmentInfo)]	
	chartPlots(y@enrichmentInfo, chartType=chartType, specifiedCols=orderby[1], ylab=orderby[1], newWindow=TRUE, ...)
}

