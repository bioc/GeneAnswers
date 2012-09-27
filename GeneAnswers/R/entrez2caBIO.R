`entrez2caBIO` <- 
function(lls) {
	lls <- unlist(lls)
	singleMap <- function (ll) {
		xmlLink <- paste('http://cabioapi.nci.nih.gov/cabio43/GetXML?query=gov.nih.nci.cabio.domain.Gene&gov.nih.nci.common.domain.DatabaseCrossReference[@sourceType=Entrez gene][@crossReferenceId=',ll, ']', sep='')
		return(.getcaBIOIDInfo(xmlLink))
	}
	print('Mapping Entrez IDs to caBIO gene IDs ...') 
	temp <- lapply(lls, singleMap)
	if (is.null(names(lls))) names(temp) <- lls
	else names(temp) <- names(lls)
	return(temp)
}
