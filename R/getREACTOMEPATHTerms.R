`getREACTOMEPATHTerms` <- 
function(pathDBIDs) {
	if (length(pathDBIDs) == 0) return(NULL)
	require(biomaRt)
	mart <- useMart("REACTOME")
	datasets <- listDatasets(mart)
	pathway<-useDataset("pathway",mart)
	pathName <- .matrix2list(as.matrix(getBM(attributes=c('pathway_db_id', '_displayname'), filters='pathway_db_id_list', values=pathDBIDs, mart=pathway)), verbose=FALSE)
	return(pathName[pathDBIDs])
}