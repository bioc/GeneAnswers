`getConnectedGraph` <- 
function(graphIDs, idType=c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'GeneInteraction', 'Customized'), edgeM=NULL, limitedLayers=FALSE, layers=6, treeMergeFilter=FALSE,
						searchAll=FALSE, showAllNodes=FALSE, directed=FALSE, direction=c('up', 'down', 'both'), filterGraphIDs=NULL, filterLayer=1, verbose=TRUE, ...) {
	idType <- match.arg(idType)
	direction <- match.arg(direction)
	#if ((length(graphIDs) < 2) & (is.null(filterGraphIDs))) stop('Input IDs should be more than 1 without filterGraphIDs! Aborting ...')
	if (limitedLayers & (layers <= 0)) stop('Specified layer can not be less than 1! Aborting ...')
	if (treeMergeFilter) {
		if (verbose) print('Search Tree Merge with filterGraphIDs!')
		if (is.vector(filterGraphIDs)) tempFilter=filterGraphIDs
		if (is.matrix(filterGraphIDs) | is.data.frame(filterGraphIDs)) tempFilter=filterGraphIDs[,1]
	} else {
		if (verbose) print('Search Tree Merge without filterGraphIDs!')
		tempFilter=NULL
	}
	
	if (directed) {
		directionVector <- switch(direction,
			'up'=TRUE,
			'down'=FALSE,
			'both'=c(TRUE, FALSE))
	} else {
		directionVector = TRUE
	}
	stopLayer <- rep(0, length(directionVector))
	connected <- rep(FALSE, length(directionVector))
	for (i in 1:length(directionVector)) {
		allLayers <- as.list(graphIDs)
		names(allLayers) <- graphIDs
		currentLayers <- graphIDs
		if (verbose) print('Searching hubs ...')
		while (TRUE) {
			if (idType %in% c('GO', 'GO.BP', 'GO.CC', 'GO.MF')) currentLayers <- getNextGOIDs(currentLayers, GOType=idType, UP=directionVector[i], filterGOIDs=tempFilter)
			else currentLayers <- getSingleLayerGraphIDs(currentLayers, edgeM, UP=directionVector[i], filterGraphIDs=tempFilter)
			currentLayers <- lapply(currentLayers, unique)
			currentLayers <- currentLayers[which(sapply(currentLayers, length) != 0)]
			currentLayers <- currentLayers[!(is.na(currentLayers))]
			if (length(currentLayers) < 1 ) break
			knownNodes <- unique(c(names(allLayers), unlist(allLayers)))
			allLayers <- lapply(allLayers , function(x, y) {return(unique(unlist(c(x, y[names(y) %in% x]))))}, currentLayers)
			if (length(allLayers) > 1) allLayers <- .treeMerge(allLayers)
			stopLayer[i] <- stopLayer[i] + 1
			currentLayers <- unique(unlist(currentLayers))
			currentLayers <- currentLayers[!(currentLayers %in% knownNodes)]
			if (length(allLayers) == 1 & (!searchAll)) break 
			if ((limitedLayers) & (stopLayer[i] >= layers)) break
		}
		if (length(allLayers) == 1) connected[i] <- TRUE
	}
	
	stopLayer = max(stopLayer)
	if (stopLayer == 0) return(print('No connection between given IDs!'))
	else {
		if (verbose) {
			if (all(!connected)) print(paste('The Graph without any nodes removal is not a connected one. Search ', stopLayer,' layer(s). Drawing network ...', sep=''))
			else print(paste('The Graph without any nodes removal is a connected one. Search ', stopLayer,' layer(s). Drawing network ...', sep=''))
		}
		if (limitedLayers) stopLayer = layers
		if (showAllNodes) filterLayer=stopLayer
		else filterLayer=filterLayer 
		return(buildNet(graphIDs, idType=idType, edgeM=edgeM, layers=stopLayer, filterGraphIDs=filterGraphIDs, verbose=verbose, filterLayer=filterLayer, directed=directed, direction=direction, ...))
	}
}
