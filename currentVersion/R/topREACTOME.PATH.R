`topREACTOME.PATH` <-
function(x, catTerm=TRUE, keepID=TRUE, ...) {
	if (length(grep('REACTOME.PATH', x@categoryType)) < 1) stop(paste('The input geneAnswers categoryType is not REACTOME.PATH but ', x@categoryType, '. stop function!'))
	if (catTerm) {
		if (keepID) rownames(x@enrichmentInfo) <- paste(getREACTOMEPATHTerms(rownames(x@enrichmentInfo)), '::', rownames(x@enrichmentInfo), sep='')
		else rownames(x@enrichmentInfo) <- getREACTOMEPATHTerms(rownames(x@enrichmentInfo))
	}
	return(topCategory(x, ...))	
}