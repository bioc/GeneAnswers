`geneAnswersReadable` <-
function(x, catTerm=TRUE, geneSymbol=TRUE, strict=FALSE, verbose=TRUE, missing=c('name', 'keep', 'remove')) {
	missing <- match.arg(missing)
	y <- x
	if (geneSymbol) {
		if (verbose) print('Converting geneInput ...')
		y@geneInput[,1] <- getSymbols(as.vector(x@geneInput[,1]), x@annLib, strict=strict, missing=missing)		
	}
 	if (verbose) print('Converting genesInCategory ...')
	if (geneSymbol) y@genesInCategory <- lapply(x@genesInCategory, getSymbols, x@annLib, strict=strict, missing=missing)
	if (catTerm) names(y@genesInCategory) <- getCategoryTerms(names(x@genesInCategory), x@categoryType, strict=strict, missing=missing)
	if (catTerm) {
		if (verbose) print('Converting enrichmentInfo rownames ...')
		rownames(y@enrichmentInfo) <- getCategoryTerms(rownames(x@enrichmentInfo), x@categoryType, strict=strict, missing=missing)
	}
	if (!is.null(y@geneExprProfile) & geneSymbol) {
		if (verbose) print('Converting geneExprProfile rownames ...')
		rownames(y@geneExprProfile) <- getSymbols(rownames(x@geneExprProfile), x@annLib, strict=strict, missing=missing)
	}
	return(y)
}

