`getDOLiteTerms` <-
function(DOLiteIDs) {
	if (is.null(DOLiteIDs)) stop('Input is NULL!')
	data('DOLiteTerm', package='GeneAnswers')
	DOLiteIDs <- as.vector(DOLiteIDs)
	temp <- DOLiteTerm[DOLiteIDs]
	names(temp) <- DOLiteIDs
	return(as.list(temp))
}

