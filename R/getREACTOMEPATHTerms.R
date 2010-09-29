`getREACTOMEPATHTerms` <- 
function(pathDBIDs, allowNA=TRUE) {
	if (length(pathDBIDs) == 0) return(NULL)
	require(biomaRt)
	if ("REACTOME" %in% toupper(listMarts()[,'biomart']))  {
		mart <- useMart("REACTOME")
		datasets <- listDatasets(mart)
		pathway<-useDataset("pathway",mart)
		returnTable <- getBM(attributes=c('pathway_db_id', '_displayname'), filters='pathway_db_id_list', values=pathDBIDs, mart=pathway)
		pathName <- .matrix2list(matrix(c(returnTable[,'pathway_db_id'], returnTable[,'_displayname']), ncol=dim(returnTable)[2], dimnames=list(c(1:dim(returnTable)[1]), colnames(returnTable))), verbose=FALSE)
		if (allowNA) {
			legalStr <- sapply(pathName, nchar, allowNA=allowNA)
			return(c(lapply(pathName[which(is.na(legalStr))], function(x) return("The name of this pathway contains special characters and is replaced!")), pathName[which(!is.na(legalStr))])[pathDBIDs]) 
		} else {
			print("If you see a warning like : input string  is invalid in this locale, that means names of pathways might contain special characters!") 
			return(pathName[pathDBIDs])
		}
		
	} else {
		stop('REACTOME service is currently not available! Please try later.') 
	}
	
}