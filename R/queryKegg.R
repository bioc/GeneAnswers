`queryKegg` <- 
function(x, ...) {
	require(KEGGREST)
	results <- list()
	step <- 10
	for (i in 1:((length(x)-1)%/%10+1)) {
		if(i == ((length(x)-1)%/%10+1)) step <- length(x)-(i-1)*10;
 		# temp <- keggGet(x[(i-1)*10+(1:step)], ...)
 		temp <- tryCatch(keggGet(x[(i-1)*10+(1:step)], ...), error=function(e) NULL)
		names(temp) <- sapply(temp, function(x) return(x$ENTRY))
		#names(temp) <- paste(organism, sapply(temp, function(x) return(x$ENTRY)), sep="")
		#names(temp) <- x[(i-1)*10+(1:step)]
		results <- c(results, temp[sapply(temp, length) != 0])
	}
	return(results)
}