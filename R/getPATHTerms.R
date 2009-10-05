`getPATHTerms` <-
function(pathIDs) {
	if (length(pathIDs) == 0) return(NULL)
	require(KEGG.db)
	pathIDs <- gsub('[^\\d]', '', pathIDs, perl=TRUE)
	pathName <- lookUp(pathIDs, 'KEGG', 'PATHID2NAME') 
	return(pathName[pathIDs])
}

