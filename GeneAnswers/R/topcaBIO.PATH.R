`topcaBIO.PATH` <-
function(x, catTerm=TRUE, keepID=TRUE, ...) {
	if (length(grep('CABIO.PATH', x@categoryType)) < 1) stop(paste('The input geneAnswers categoryType is not CABIO.PATH but ', x@categoryType, '. stop function!'))
	if (catTerm) {
		if (keepID) rownames(x@enrichmentInfo) <- paste(getcaBIOPATHTerms(rownames(x@enrichmentInfo)), '::', rownames(x@enrichmentInfo), sep='')
		else rownames(x@enrichmentInfo) <- getcaBIOPATHTerms(rownames(x@enrichmentInfo))
	}
	return(topCategory(x, ...))	
}