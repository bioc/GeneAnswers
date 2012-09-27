`geneAnswersSort` <-
function(x, sortBy = c('geneNum', 'pvalue', 'foldChange', 'oddsRatio', 'correctedPvalue')) {
	sortBy <- match.arg(sortBy)
	if ((sortBy == 'correctedPvalue') & !('fdr p value' %in% colnames(x@enrichmentInfo))) stop('input geneAnswer class does not contain corrected p value!!!')
	orderby <- switch(sortBy,
		'geneNum'= c('genes in Category', 'TRUE'),
		'pvalue' = c('p value', 'FALSE'),
		'foldChange' = c('fold of overrepresents', 'TRUE'),
		'oddsRatio' = c('odds ratio', 'FALSE'), 
		'correctedPvalue' = c('fdr p value', 'FALSE'))
	y <- x
	y@enrichmentInfo <-	x@enrichmentInfo[order(x@enrichmentInfo[, orderby[1]], decreasing=as.logical(orderby[2])), ]
	y@genesInCategory <- x@genesInCategory[c(rownames(y@enrichmentInfo), names(x@genesInCategory)[!(names(x@genesInCategory) %in% rownames(x@enrichmentInfo))])]
	return(y)
}

