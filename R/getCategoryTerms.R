`getCategoryTerms` <-
function(catIDs, catType, strict=FALSE, missing=c('name', 'keep', 'remove'), nameLength='all', addID=FALSE) {
	missing=match.arg(missing)
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
		if (tolower(nameLength) != 'all') {
			if ((length(nameLength) == 1) & is.numeric(nameLength)) {
				tempNames <- names(temp)
				dots <- rep('...', length(temp))
				dots[which(nchar(temp) <= nameLength)] <- ''
				temp <- paste(sapply(temp, substr, 1, nameLength), dots, sep='')
				names(temp)  <- tempNames
			} else {
				stop('nameLength can not be recognized!')
			}
		}
		if (addID) {
			tempNames <- names(temp)
			temp <- paste(temp, names(temp), sep='::')
			names(temp) <- tempNames
		}
		return(temp)
	} else {
		print('Slot categoryType is not recognized! No mapping ...')
		return(catIDs)
	}
}

