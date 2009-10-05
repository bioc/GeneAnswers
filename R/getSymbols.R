`getSymbols` <-
function (geneIDs, data, strict=FALSE, missing=c('name', 'keep', 'remove')) {
	missing <- match.arg(missing)
	temp <- getSYMBOL(geneIDs, data)
	if (NA %in% temp) {
		print('Warning: some genes do not have valid symbols!')
		if (strict) stop('Interrupt conversion!')
	}
	if (missing != 'keep') {
		if (missing == 'name') temp[temp %in% NA] <- names(temp[temp %in% NA])
		else temp <- temp[!(temp %in% NA)]
	}
	filteredGeneIDs <- geneIDs[geneIDs %in% names(temp)]
	return(temp[filteredGeneIDs])
}

