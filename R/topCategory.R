`topCategory` <-
function(inputX, orderby=c('geneNum', 'pvalue', 'foldChange', 'oddsRatio', 'correctedPvalue'), top=5, file=FALSE, fileName='topCategory.txt') {
	orderby <- match.arg(orderby)
	if (is.character(top) & toupper(top) != 'ALL') stop('top can not be recognized!')
	if (is.numeric(top)) top <- min(dim(inputX@enrichmentInfo)[1], top)
	else {
		if (toupper(top) == 'ALL') top <- min(dim(inputX@enrichmentInfo)[1], 30)
		print(paste(top, 'categories will be displayed! All categories can be saved in file'))
	}
	if ((orderby == 'correctedPvalue') & !('fdr p value' %in% colnames(inputX@enrichmentInfo))) stop('input geneAnswer class does not contain fdr p value!!!')
	orderby <- switch(orderby,
		'geneNum'= c('genes in Category', 'TRUE'),
		'pvalue' = c('p value', 'FALSE'),
		'foldChange' = c('fold of overrepresents', 'TRUE'),
		'oddsRatio' = c('odds ratio', 'FALSE'), 
		'correctedPvalue' = c('fdr p value', 'FALSE'))
	inputX@enrichmentInfo <- inputX@enrichmentInfo[order(inputX@enrichmentInfo[, orderby[1]], decreasing=as.logical(orderby[2])), ]
	displayCols <- c('genes in Category', 'p value')
	if ('fdr p value' %in% colnames(inputX@enrichmentInfo)) displayCols <- c(displayCols, 'fdr p value')
	if (file) {
		outputM <- cbind(as.data.frame(matrix(rownames(inputX@enrichmentInfo[1:top,]))), inputX@enrichmentInfo[1:top, displayCols])
		names(outputM) <- c('Category', displayCols)
		write.table(outputM, file=fileName, sep='\t', row.names=F, col.names=T, append=F)
		print(paste('File', fileName, 'is successfully generated!'))
	}
	else print(inputX@enrichmentInfo[1:top, displayCols])
}

