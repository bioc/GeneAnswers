`getCategoryTerms` <-
function(catIDs, catType, strict=FALSE, missing=c('name', 'keep', 'remove')) {
	if (toupper(catType) %in% c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'DOLITE', 'KEGG')) {
		catTerms = switch(toupper(catType), 
		'GO'=.getGOTerms(catIDs),
		'GO.BP'=.getGOTerms(catIDs),
		'GO.CC'=.getGOTerms(catIDs),
		'GO.MF'=.getGOTerms(catIDs),
		'DOLITE'=getDOLiteTerms(catIDs),
		'KEGG'=getPATHTerms(catIDs))
		if (NA %in% unlist(catTerms)) {
			print('Warning: some category IDs do not have names!')
			if (strict) stop('Interrupt conversion!')
		}
		temp <- unlist(catTerms[catIDs])
		if (missing != 'keep') {
			if (missing == 'name') temp[temp %in% NA] <- names(temp[temp %in% NA])
			else temp <- temp[!(temp %in% NA)]
		}
		return(temp)
	} else {
		print('Slot categoryType is not recognized! No mapping ...')
		return(catIDs)
	}
}

