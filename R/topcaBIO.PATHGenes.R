`topcaBIO.PATHGenes` <-
function(x, catTerm=TRUE, keepID=TRUE, geneSymbol=TRUE, ...) {
	stop("Due to termination of caBig, this function is removed in this version!")
#	if (length(grep('CABIO.PATH', x@categoryType)) < 1) stop(paste('The input geneAnswers categoryType is not CABIO.PATH but ', x@categoryType, '. stop function!'))
#	if (catTerm) {
#		if (keepID) {
#			rownames(x@enrichmentInfo) <- paste(unlist(getcaBIOPATHTerms(rownames(x@enrichmentInfo))), '::', rownames(x@enrichmentInfo), sep='') 
#			names(x@genesInCategory) <- paste(unlist(getcaBIOPATHTerms(names(x@genesInCategory))), '::', names(x@genesInCategory), sep='')
#		} else {
#			rownames(x@enrichmentInfo) <- getcaBIOPATHTerms(rownames(x@enrichmentInfo))
#			names(x@genesInCategory) <- getcaBIOPATHTerms(names(x@genesInCategory))
#		}
#	} 
#	if (geneSymbol) {
#		if (is.null(x@annLib)) stop('annotation library is not specified! original gene IDs will be kept.')
#		else {
#			x@genesInCategory <- lapply(x@genesInCategory, getSymbols, x@annLib)
#			x@geneInput[,1] <- getSymbols(x@geneInput[,1], x@annLib)
#		}
#	}
#	return(topCategoryGenes(x, ...))
}

