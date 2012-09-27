`getcaBIOPATHTerms` <- 
function(caBIOPATHIDs) {
	.getPathTerms <- function(x) {
		if (is.na(x)) return(NA)
		else return(paste(x[['displayValue']], '@', x[['source']], sep=''))
	}
	op <- options()
	options(warn=-1)
	result <- lapply(.getcaBIOPATHInfo(caBIOPATHIDs), .getPathTerms)[caBIOPATHIDs]
	options(op)
	return(result)
}