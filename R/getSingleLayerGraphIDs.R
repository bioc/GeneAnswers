`getSingleLayerGraphIDs` <- 
function(graphIDs, edgeM, remove=TRUE, filterGraphIDs=NULL, directed=FALSE, UP=TRUE) { 
	singleLayerGraphIDs <- lapply(graphIDs, .getSingleGraphIDs, edgeM, filterGraphIDs=filterGraphIDs, directed=directed, UP=UP)
	names(singleLayerGraphIDs) <- as.character(graphIDs)
	if (remove) singleLayerGraphIDs <- singleLayerGraphIDs[which(sapply(singleLayerGraphIDs, length) != 0)]
	if (length(singleLayerGraphIDs) > 0) return(singleLayerGraphIDs)
	else return(NULL)
}