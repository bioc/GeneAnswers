`getPATHTerms` <-
function(pathIDs) {
	if (length(pathIDs) == 0) return(NULL)
	#require(KEGG.db)
	##pathIDs <- gsub('[^\\d]', '', pathIDs, perl=TRUE)
	#pathName <- lookUp(pathIDs, 'KEGG', 'PATHID2NAME') 
	pathName <- lapply(queryKegg(pathIDs), function(x) return(x$NAME))
	return(pathName[pathIDs])
}
