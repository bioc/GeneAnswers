`getNextGOIDs` <- 
function (GOIDs, GOType=c('GO', 'GO.BP', 'GO.CC', 'GO.MF'), remove=TRUE, filterGOIDs=NULL, UP=TRUE) {
	if (length(GOIDs) == 0) return(NULL)
	else {
		GOType <- match.arg(GOType)
		require(GO.db)
		if (UP) direction <- 'PARENTS'
		else direction <- 'CHILDREN'
		GOCat <- gsub("^GO.", '', GOType)
		if ("GO" %in% GOCat) GOCat <- c('BP', 'CC', 'MF')
		if (!('GO' %in% GOType)) GOIDs <- GOIDs[GOIDs %in% aqListGOIDs(GOCat)]
		if (length(GOIDs) == 0) return(NULL)
		parentsGOIDs <- list()
		if (('GO' %in% GOType) | ('GO.BP' %in% GOType))	parentsGOIDs <- c(parentsGOIDs, lookUp(GOIDs, "GO", paste("BP", direction, sep='')))
		if (('GO' %in% GOType) | ('GO.CC' %in% GOType))	parentsGOIDs <- c(parentsGOIDs, lookUp(GOIDs, "GO", paste("CC", direction, sep='')))
		if (('GO' %in% GOType) | ('GO.MF' %in% GOType))	parentsGOIDs <- c(parentsGOIDs, lookUp(GOIDs, "GO", paste("MF", direction, sep='')))
		parentsGOIDs <- lapply(parentsGOIDs, unique)
		if (!is.null(filterGOIDs)) parentsGOIDs <- lapply(parentsGOIDs, intersect, filterGOIDs)
		if (remove) {
			parentsGOIDs <- parentsGOIDs[which(sapply(parentsGOIDs, length) != 0)]
			parentsGOIDs <- parentsGOIDs[!(is.na(parentsGOIDs))]
		}
		if (length(parentsGOIDs) > 0) return(parentsGOIDs)
		else return(NULL)
	}
}