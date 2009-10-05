`getDOLiteTerms` <-
function(DOLiteIDs) {
	if (is.null(DOLiteIDs)) stop('Input is NULL!')
	return(as.list(DOLiteTerm[DOLiteIDs]))
}

