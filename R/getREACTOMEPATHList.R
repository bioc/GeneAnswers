`getREACTOMEPATHList` <- 
function(geneVector, lib) {
	libname <- sub('\\.db', '', lib)
	idType <- switch(sub('org.*[:.:]', '', libname), 'eg'='EG', 'tair'='TAIR', 'ORF')
	require(biomaRt)
	mart <- useMart("REACTOME")
	datasets <- listDatasets(mart)
	pathway<-useDataset("pathway",mart)
	if (idType == 'EG') {
		returnType <- c("pathway_db_id", "referencedatabase_entrez_gene")
		queryType <- 'referencednasequence_entrez_id_list'
	} else {
		if (idType == 'TAIR') {
			returnType <- c("pathway_db_id", "referencedatabase_tigr")
			queryType <- 'referencepeptidesequence_tigr_id_list'
		}
		else stop('ORF is not currently supported by biomaRt!')
	}
	print('Start to query biomaRt to retrieve related pathways ...')
	paths <- unique(as.matrix(getBM(attributes=returnType, filters=queryType, values=geneVector, mart=pathway))[,'pathway_db_id'])
	print('Start to query biomaRt to retrieve all genes for the related pathways ...')
	temp <- getBM(attributes=returnType, filters='pathway_db_id_list', values=paths, mart=pathway)
	return(lapply(.matrix2list(as.matrix(temp[!is.na(temp[,returnType[2]]),]), verbose=FALSE), as.character))
}