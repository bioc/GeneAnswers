`getMultiLayerGraphIDs` <- 
function(graphIDs, idType=c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'GeneInteraction', 'Customized'), edgeM=NULL, layers=1, 
								filterGraphIDs=NULL, filterLayer=0, UP=TRUE, directed=FALSE) {
	if (!is.numeric(layers)) stop('Layers should be an integer!')
	if ((filterLayer > layers) | (!(is.null(filterGraphIDs)) & (filterLayer < 1))) stop('filter layer should not be more than the maximal layers or less than 1!')
	if (is.null(filterGraphIDs) & (filterLayer > 0)) stop('No filter graphIDs! Can not set filter layer ...')
	if (!(is.null(filterGraphIDs))) filterGraphIDs <- as.character(filterGraphIDs)
	idType <- match.arg(idType)
	multiLayerGraphIDs <- c()
	tempGraphIDs <- as.character(graphIDs)
	layersLength <- c()
	layers <- layers + 1
	for (i in 1:layers) {
		if (i == layers) tempFilter <- unique(c(names(multiLayerGraphIDs), unlist(multiLayerGraphIDs)))
		else {
			if (is.null(filterGraphIDs)) {
				tempFilter <- NULL
				filterLayer <- 0
			} else {
				if (i < filterLayer) tempFilter <- NULL
				else tempFilter <- unique(c(filterGraphIDs, names(multiLayerGraphIDs), unlist(multiLayerGraphIDs)))
			}
		}
		
		if (idType %in% c('GO', 'GO.BP', 'GO.CC', 'GO.MF')) temp <- getNextGOIDs(tempGraphIDs, GOType=idType, filterGOIDs=tempFilter, UP=UP)
		else {
			if (idType %in% c('GeneInteraction', 'Customized')) temp <- getSingleLayerGraphIDs(tempGraphIDs, edgeM, filterGraphIDs=tempFilter, UP=UP)
			else stop('The given idType is not supported!')
		}
		
		if ((length(layersLength) > 0) & !directed & (i <=  filterLayer)) {
			l <- layersLength[length(layersLength)]
			tempL <- multiLayerGraphIDs[(length(multiLayerGraphIDs)-l+1):length(multiLayerGraphIDs)]
			tempLM <- .list2matrix(tempL)
			tempCM <- .list2matrix(temp)
			tempDLM <- rbind(tempLM, cbind(tempCM[,2], tempCM[,1]))
			tempCM <- tempCM[!duplicated(tempDLM)[-(1:dim(tempLM)[1])], ]
			if (!is.matrix(tempCM)) tempCM <- matrix(tempCM, ncol=2)
			temp <- .matrix2list(tempCM)
		}
		
		if ((i == filterLayer) & (length(tempGraphIDs) > length(temp))) {
		    j <- filterLayer - 1
			tempList <- temp
			tempMultiLayerGraphIDs <- c()
			while (j > 0) {
				tempParentsLayers <- multiLayerGraphIDs[(length(multiLayerGraphIDs) - layersLength[j] + 1):length(multiLayerGraphIDs)]
				tempParentsLayersNodes <- unlist(tempParentsLayers)
				tempList <- lapply(tempParentsLayers, intersect, unique(c(names(tempList), names(multiLayerGraphIDs), unlist(tempList), filterGraphIDs, tempParentsLayersNodes[duplicated(tempParentsLayersNodes)])))
				tempList <- tempList[which(sapply(tempList, length) != 0)]
				layersLength[j] <- length(tempList)
				tempMultiLayerGraphIDs <- c(tempList, tempMultiLayerGraphIDs)
				j <- j - 1
				multiLayerGraphIDs <- multiLayerGraphIDs[!(names(multiLayerGraphIDs) %in% names(tempParentsLayers))]
			}
			multiLayerGraphIDs <- tempMultiLayerGraphIDs
		}                                                                              
		if (is.null(temp)) break
		
		tempGraphIDs <- unique(unlist(temp))
		tempGraphIDs <- tempGraphIDs[!(tempGraphIDs %in% unique(c(names(multiLayerGraphIDs), unlist(multiLayerGraphIDs), names(temp))))]
		multiLayerGraphIDs <- c(multiLayerGraphIDs, temp)
		
		layersLength <- c(layersLength, length(temp))
		#parentsGenes <- names(multiLayerGraphIDs)
		
		#tempGraphIDs <- tempGraphIDs[!(tempGraphIDs %in% parentsGenes)]
	}
	if (i < layers) selfClosure <- TRUE
	else selfClosure <- FALSE
	return(c(list(selfClosure), list(layersLength=layersLength), multiLayerGraphIDs))
}
