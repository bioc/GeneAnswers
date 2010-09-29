`topCategoryGenes` <-
function(inputX, orderby=c('geneNum', 'pvalue', 'foldChange', 'oddsRatio', 'correctedPvalue'), top=5, genesOrderBy = 1, decreasing=FALSE, topGenes = 5, file=FALSE, fileName='topCategoryGenes.txt') {
	orderby <- match.arg(orderby)
	allCats <- FALSE
	allGenes <- FALSE
	if (is.character(top) & toupper(top) != 'ALL') stop('top can not be recognized!')
	if (is.numeric(top)) top <- min(dim(inputX@enrichmentInfo)[1], top)
	else if (toupper(top) == 'ALL') {
		allCats <- TRUE
		top <- min(dim(inputX@enrichmentInfo)[1], 6)
		print(paste(top, 'categories will be displayed! All categories are saved in file'))
	}
	if ((orderby == 'correctedPvalue') & !('fdr p value' %in% colnames(inputX@enrichmentInfo))) stop('input geneAnswer class does not contain fdr p value!!!')
	orderby <- switch(orderby,
		'geneNum'= c('genes in Category', 'TRUE'),
		'pvalue' = c('p value', 'FALSE'),
		'foldChange' = c('fold of overrepresents', 'TRUE'),
		'oddsRatio' = c('odds ratio', 'FALSE'), 
		'correctedPvalue' = c('fdr p value', 'FALSE'))
	inputX@enrichmentInfo <- inputX@enrichmentInfo[order(inputX@enrichmentInfo[, orderby[1]], decreasing=as.logical(orderby[2])), ]
	if (is.character(topGenes) & toupper(topGenes) != 'ALL') stop('topGenes can not be recognized!')
	if (toupper(genesOrderBy) == 'NONE') sortGenes <- FALSE
	else sortGenes <- TRUE
	if (allCats) top <- dim(inputX@enrichmentInfo)[1]
	if (file) writeLines(paste(paste(top, ' categories genes table', sep=''), '', paste(c('Category', orderby[1], colnames(inputX@geneInput)), collapse='\t'), sep='\n'), con=fileName)
	for (i in 1:top) {
		tempTopGenes <- topGenes
		temp <-inputX@geneInput[inputX@geneInput[,1] %in% inputX@genesInCategory[[rownames(inputX@enrichmentInfo)[i]]],]
		if (sortGenes) {
			if (dim(inputX@geneInput)[2] == 1) {
				temp <- as.data.frame(sort(temp), stringsAsFactors=FALSE)
				colnames(temp) <- colnames(inputX@geneInput)
			} 
			else temp <- temp[order(temp[, genesOrderBy], decreasing=decreasing),]
		}
			
		if (is.numeric(tempTopGenes)) tempTopGenes <- min(dim(temp)[1], tempTopGenes)
		else {
			if (toupper(tempTopGenes) == 'ALL') {
				tempTopGenes <- dim(temp)[1]
				allGenes <- TRUE
			}
		}
		if (file) {
			writeDataFrame <- cbind(as.data.frame(matrix(c(rep(rownames(inputX@enrichmentInfo)[i], length=tempTopGenes), rep(inputX@enrichmentInfo[i, orderby[1]], length=tempTopGenes)), nrow = tempTopGenes, ncol=2, byrow=FALSE)), temp[1:tempTopGenes, ])
			names(writeDataFrame) <- c('Category', orderby[1], names(temp))
			write.table(writeDataFrame, file=fileName, sep='\t', row.names=F, col.names=F, append=T)
		}
		if ((allCats & (i <= 6)) | (!allCats)) {
			print(paste('******** ', rownames(inputX@enrichmentInfo)[i], '   ', orderby[1], ' : ', signif(as.numeric(inputX@enrichmentInfo[i, orderby[1]]), digits=4), ' ********', sep='')) 
			if (allGenes) print(temp[1:min(5, dim(temp)[1]), ])
			else print(temp[1:min(tempTopGenes, dim(temp)[1]), ])
		} 
	}
	if (file) print(paste('File', fileName, 'is successfully generated!'))
}

