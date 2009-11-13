`geneAnswersConceptNet` <-
function(x, colorValueColumn=NULL, centroidSize=c('pvalue', 'geneNum', 'foldChange', 'oddsRatio', 'correctedPvalue'), output=c('fixed','interactive'), 
			showCats=c(1:5), geneLayer=1, edgeM=NULL, catTerm=FALSE, geneSymbol=FALSE, catID=FALSE, nameLength='all', ...) {
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
	if (is.numeric(showCats)) {
		if (!(all(showCats %in% c(1:dim(x@enrichmentInfo)[1])))) print('Some specified categories might not be statistical significant! Only show significant categories.')
		showCats <- intersect(showCats, c(1:dim(x@enrichmentInfo)[1])) 
	} else {
		if (is.character(showCats)) {
			showCats <- intersect(showCats, rownames(x@enrichmentInfo))
			if (length(showCats) < 1) stop('specified categories can not be recognized!')
		} else stop('specified categories can not be recognized!')
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
  			names(inputList) <- getCategoryTerms(names(inputList), x@categoryType, missing='name', nameLength=nameLength, addID=catID)
  			names(scaledTemp) <- getCategoryTerms(names(scaledTemp), x@categoryType, missing='name', nameLength=nameLength, addID=catID)
		} else {
			print('Slot categoryType is not recognized! No mapping ...')
		}
	} 
	
	extendIA <- function(IA) {
		if (dim(IA)[2] > 2) tempM <- cbind(IA[,2], IA[,1], IA[,3:dim(IA)[2]])
		else tempM <- cbind(IA[,2], IA[,1])
		colnames(tempM) <- colnames(IA)
		tempM <- rbind(IA, tempM)[,1:2]
		return(tempM[!(duplicated(tempM)),])
	}
	
	if (geneLayer > 1) {
		if (is.null(x@annLib)) {
			if (is.null(edgeM)) stop("Customized database is not available!")
		} else {
			switch(x@annLib,
				'org.Hs.eg.db'=data('HsIALite', package='GeneAnswers'),
				'org.Mm.eg.db'=data('MmIALite', package='GeneAnswers'),
				'org.Rn.eg.db'=data('RnIALite', package='GeneAnswers'),
				'org.Dm.eg.db'=data('DmIALite', package='GeneAnswers'))
			edgeM <- switch(x@annLib,
				'org.Hs.eg.db'=extendIA(HsIALite),
		   		'org.Mm.eg.db'=extendIA(MmIALite),
				'org.Rn.eg.db'=extendIA(RnIALite),
				'org.Dm.eg.db'=extendIA(DmIALite))
		}
		IAgenes <- getMultiLayerGraphIDs(unique(unlist(inputList)), idType='GeneInteraction', edgeM=edgeM, layers=(geneLayer-1), filterGraphIDs=x@geneInput[,1], filterLayer=1)
		IAgenes <- IAgenes[-1:-2]
		if (geneSymbol & (length(IAgenes) > 1)) {
			IAgenes <- lapply(IAgenes, getSymbols, x@annLib, missing='name')
			names(IAgenes) <- getSymbols(names(IAgenes), x@annLib, missing='name')
		}
	}
	
	if (geneSymbol) {
		if (!is.null(inputXValue)) names(inputXValue) <- getSymbols(names(inputXValue), x@annLib, missing='name')
		inputList <- lapply(inputList, getSymbols, x@annLib, missing='name')
	}
	
	if (geneLayer > 1) {
		if (length(IAgenes) > 1) inputList <- c(inputList, IAgenes)
	} 
	
	geneConceptNet(inputList, lengthOfRoots=length(scaledTemp), inputValue = inputXValue[unique(c(unlist(inputList),names(inputList)[-(1:length(scaledTemp))]))], centroidSize=scaledTemp, output=output, ...)
	return(invisible(x))
}

