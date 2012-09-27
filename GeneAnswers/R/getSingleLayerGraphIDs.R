`getSingleLayerGraphIDs` <- 
function(graphIDs, edgeM, remove=TRUE, filterGraphIDs=NULL, UP=TRUE) { 
	singleLayerGraphIDs <- lapply(graphIDs, .getSingleGraphIDs, edgeM, filterGraphIDs=filterGraphIDs, UP=UP)
	names(singleLayerGraphIDs) <- as.character(graphIDs)
	if (remove) singleLayerGraphIDs <- singleLayerGraphIDs[which(sapply(singleLayerGraphIDs, length) != 0)]
	if (length(singleLayerGraphIDs) > 0) return(singleLayerGraphIDs)
	else return(NULL)
}