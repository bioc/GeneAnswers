`topcaBIO.PATH` <-
function(x, catTerm=TRUE, keepID=TRUE, ...) {
	stop("Due to termination of caBig, this function is removed in this version!")
	#if (length(grep('CABIO.PATH', x@categoryType)) < 1) stop(paste('The input geneAnswers categoryType is not CABIO.PATH but ', x@categoryType, '. stop function!'))
	#if (catTerm) {
	#	if (keepID) rownames(x@enrichmentInfo) <- paste(getcaBIOPATHTerms(rownames(x@enrichmentInfo)), '::', rownames(x@enrichmentInfo), sep='')
	#	else rownames(x@enrichmentInfo) <- getcaBIOPATHTerms(rownames(x@enrichmentInfo))
	#}
	#return(topCategory(x, ...))	
}