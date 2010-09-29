`getConceptTable` <-
function (gAList, topCat=10, items=c('both', 'geneNum', 'pvalue'), sortBy = c('pvalue', 'geneNum', 'foldChange', 'oddsRatio', 'correctedPvalue'), catTerm=TRUE, strict=FALSE) {
	items <- match.arg(items)
	sortBy <- match.arg(sortBy)
	
	if (is.numeric(topCat) & length(topCat) == 1) {
		topCatFun <- function (x, top=topCat) { top <- min(top, dim(x@enrichmentInfo)[1]); return(rownames(geneAnswersSort(x, sortBy=sortBy)@enrichmentInfo)[1:top]) }
		catList <- lapply(gAList, topCatFun, top=topCat)
		categories <- unique(unlist(catList))
	} else {
		if (is.character(topCat)) {
			categories <- topCat
		} else {
			pickCat <- function (x, pick=topCat) { pick <- intersect(pick, c(1:dim(x@enrichmentInfo)[1])); return(rownames(geneAnswersSort(x, sortBy=sortBy)@enrichmentInfo)[pick]) }
			catList <- lapply(gAList, pickCat, pick=topCat)
			categories <- unique(unlist(catList))
		}
	}
	
    if (NA %in% categories) stop('Some given categories can not be found in GeneAnswers instances!!!')

	if (items == 'both') {
		catTable <- matrix('0 (1)', nrow=length(categories), ncol=length(gAList))
  		rownames(catTable) <- categories
  		colnames(catTable) <- names(gAList) 
  		for (i in 1:length(gAList)) {
  			temp <- unlist(lapply(getGenesInCategory(gAList[[i]]), length))
  			temp <- temp[names(temp) %in% categories]
  			catTable[names(temp),i] <- paste(temp, ' (p > ',gAList[[i]]@pvalueT, ')', sep='') 
  		}
  	}
  	else {
  		if (items == 'geneNum') catTable <- matrix(0, nrow=length(categories), ncol=length(gAList))
  		else catTable <- matrix(1, nrow=length(categories), ncol=length(gAList))
  		rownames(catTable) <- categories
  		colnames(catTable) <- names(gAList) 
  		if (items == 'geneNum') {
  			for (i in 1:length(gAList)) {
  				temp <- unlist(lapply(getGenesInCategory(gAList[[i]]), length))
  				temp <- temp[names(temp) %in% categories]
  				catTable[names(temp),i] <- temp
  			}
  		}
  	}

  	indexTable <- matrix(1, nrow=length(categories), ncol=length(gAList))
  	rownames(indexTable) <- categories
  	colnames(indexTable) <- names(gAList)

  	for (i in 1:length(gAList)) {
  		temp <- gAList[[i]]@enrichmentInfo[rownames(gAList[[i]]@enrichmentInfo) %in% categories, ]
  		if (items == 'both') tempDF <- paste(temp[,'genes in Category'], ' (', signif(temp[, 'p value'], digits=3), ')', sep='')
  		else tempDF <- temp[,'genes in Category']
  		names(tempDF) <- rownames(temp)
  		catTable[rownames(temp),i] <- tempDF
  		indexTable[rownames(temp), i] <- temp[, 'p value']
  	}
  	catTable <- rbind(catTable, unlist(lapply(gAList, function(x) return(length(getGeneInput(x)[,1])))))
    if (catTerm) {
		rownames(catTable)[1:(dim(catTable)[1]-1)] <- getCategoryTerms(rownames(catTable)[1:(dim(catTable)[1]-1)], getCategoryType(gAList[[1]]), strict=strict, missing='name')
		rownames(indexTable) <- rownames(catTable)[1:(dim(catTable)[1]-1)]
	} 
  	rownames(catTable)[length(rownames(catTable))] <- 'Genes / Group'
  	result <- list(as.data.frame(catTable, stringsAsFactors =FALSE), indexTable)
  	names(result) <- c('CategoriesTable', 'IndexTable')
  	return(result)
}