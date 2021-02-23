#`getPATHList` <-
#function(geneVector, lib) {
#	return(getCategoryList(geneVector, lib, 'PATH'))
#}

`getPATHList` <-
function(geneVector, lib) {
	organism <- .mappingOrganisms(lib)
	if (!is.null(organism)) {
		query <- queryKegg(paste(organism, ":", geneVector, sep=""))
		#query <- queryKegg(as.character(geneVector))
		names(query) <- NULL
		pathTerms <- unique(names(unlist(lapply(query, function(x) return(x$PATHWAY)))))
		if (!is.null(pathTerms)) {
			pathList <- sapply(queryKegg(pathTerms), function(x) return(x$GENE[seq(1,length(x$GENE)-1,length.out=length(x$GENE)/2)]))
			pathwayList <- pathList[sapply(pathList, length) != 0]
			if (!is.null(pathwayList)) {
				return(pathwayList)
			} else {
				stop('input genes do not belong to any category!!!')
			}
			return()
		}
		else {
			stop('input genes do not belong to any category!!!')
		}
	} else {
		stop('the specified annotation library can not be found!!!')
	}
}