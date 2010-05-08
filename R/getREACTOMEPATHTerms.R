`getREACTOMEPATHTerms` <- 
function(pathDBIDs) {
	if (length(pathDBIDs) == 0) return(NULL)
	require(biomaRt)
	mart <- useMart("REACTOME")
	datasets <- listDatasets(mart)
	pathway<-useDataset("pathway",mart)
	returnTable <- getBM(attributes=c('pathway_db_id', '_displayname'), filters='pathway_db_id_list', values=pathDBIDs, mart=pathway)
	pathName <- .matrix2list(matrix(c(returnTable[,'pathway_db_id'], returnTable[,'_displayname']), ncol=dim(returnTable)[2], dimnames=list(c(1:dim(returnTable)[1]), colnames(returnTable))), verbose=FALSE)
	return(pathName[pathDBIDs])
}