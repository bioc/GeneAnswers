`caBIO2entrez` <- 
function(caBIOIds) {
	caBIOIds <- unlist(caBIOIds)
	singleMap <- function(caBIOId) {
		xmlLink <- paste('http://cabioapi.nci.nih.gov/cabio43/GetXML?query=gov.nih.nci.common.domain.DatabaseCrossReference&gov.nih.nci.cabio.domain.Gene[@id=',caBIOId, ']', sep='')
		return(.getcaBIOIDInfo(xmlLink, IDinfo='direct'))
	}
	print('Mapping caBIO gene IDs to Entrez IDs ...')
	temp <- lapply(caBIOIds, singleMap)
	if (is.null(names(caBIOIds))) names(temp) <- caBIOIds
	else names(temp) <- names(caBIOIds)
	return(temp)	
}