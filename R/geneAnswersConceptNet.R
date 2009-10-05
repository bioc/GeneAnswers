`geneAnswersConceptNet` <-
function(x, colorValueColumn=NULL, centroidSize=c('geneNum', 'pvalue', 'foldChange', 'oddsRatio', 'correctedPvalue'), output=c('fixed','interactive'), showCats=c(1:5), catTerm=FALSE, geneSymbol=FALSE) {
	centroidSize <- match.arg(centroidSize)
	if ((centroidSize == 'correctedPvalue') & !('fdr p value' %in% colnames(x@enrichmentInfo))) stop('input geneAnswer class does not contain fdr p value!!!')
	#x <- geneAnswersReadable(x, catTerm=catTerm, geneSymbol=geneSymbol)
	orderby <- switch(centroidSize,
		'geneNum'= c('genes in Category', 'TRUE', 'Normal'),
		'pvalue' = c('p value', 'FALSE', '-Log10'),
		'foldChange' = c('fold of overrepresents', 'TRUE', 'Normal'),
		'oddsRatio' = c('odds ratio', 'FALSE', '-Log'), 
		'correctedPvalue' = c('fdr p value', 'FALSE', '-Log10'))
	x@enrichmentInfo <- x@enrichmentInfo[order(x@enrichmentInfo[, orderby[1]], decreasing=as.logical(orderby[2])), ]
	centroidSize <- x@enrichmentInfo[, orderby[1]]	
	names(centroidSize) <- rownames(x@enrichmentInfo)
	#if (is.character(top) & toupper(top) != 'ALL') stop('top can not be recognized!')
	#if (toupper(top) == 'ALL') top <- dim(x@enrichmentInfo)[1]
	#if (is.numeric(top)) top <- min(dim(x@enrichmentInfo)[1], top)
	if ((dim(x@geneInput)[2] == 1) | is.null(colorValueColumn)) inputXValue <- NULL
	else {
		if (is.numeric(as.numeric(try(x@geneInput[, colorValueColumn])))) {
			if (NA %in% as.numeric(x@geneInput[, colorValueColumn])) {
				print(paste('Specified ', colorValueColumn, ' does not contain valid values, No value will be assigned!'))
				inputXValue <- NULL 
			} else {
				inputXValue <- as.numeric(x@geneInput[,colorValueColumn])
				names(inputXValue) <- x@geneInput[,1]
			}
		}
	}
	if (is.numeric(showCats)) showCats <- intersect(showCats, c(1:dim(x@enrichmentInfo)[1]))
    else {
		if (is.character(showCats)) showCats <- intersect(showCats, rownames(x@enrichmentInfo))
		else stop('specified categories can not be recognized!')
	}
    inputList <- x@genesInCategory[names(x@genesInCategory) %in% names(centroidSize[showCats])]
#	newInput <- lapply(inputList, getSYMBOL, x@annLib)

#	names(newInput) <- getCategoryTerms(names(inputList), inputX@categoryType)
#	unlist(getGOTerm(names(inputList)))

#	if (!is.null(inputXValue)) names(inputXValue) <- getSYMBOL(names(inputXValue), inputX@annLib)
	temp <- centroidSize[showCats][names(inputList)]

#	names(temp) <- getCategoryTerms(names(temp), inputX@categoryType)
#	unlist(getGOTerm(names(temp)))

	scaledTemp <- switch(orderby[3],
		'Normal'=temp,
		'-Log'=-log(temp),
		'-Log10'= -log10(temp))
	if (catTerm) {
		if (x@categoryType %in% c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'DOLite', 'KEGG')) {
			names(inputList) <- getCategoryTerms(names(inputList), x@categoryType)
			names(scaledTemp) <- getCategoryTerms(names(scaledTemp), x@categoryType)
		} else {
			print('Slot categoryType is not recognized! No mapping ...')
		}
	} 
	
	if (geneSymbol) {
		if (!is.null(inputXValue)) names(inputXValue) <- lookUp(names(inputXValue), x@annLib, 'SYMBOL')
		inputList <- lapply(inputList, lookUp, x@annLib, 'SYMBOL')
	}
	geneConceptNet(inputList, inputValue = inputXValue, centroidSize=scaledTemp, output=output)
	return(invisible(x))
}

