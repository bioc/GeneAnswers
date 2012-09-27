`getSymbols` <-
function (geneIDs, data, strict=FALSE, missing=c('name', 'keep', 'remove')) {
	libname <- sub('\\.db', '', data)
	idType <- switch(sub('org.*[:.:]', '', libname), 'eg'='EG', 'tair'='TAIR', 'ORF')
	missing <- match.arg(missing)
	require(data, character.only=TRUE)
	temp <- getSYMBOL(geneIDs, data)
	if (idType == 'TAIR') {
		temp1 <- lapply(geneIDs, grep, names(temp))
		names(temp1) <- geneIDs
		temp <- temp[unlist(lapply(temp1, min))[geneIDs]]
		names(temp) <- geneIDs
	}
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

