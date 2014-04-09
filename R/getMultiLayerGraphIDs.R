`getMultiLayerGraphIDs` <- 
function(graphIDs, idType=c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'GeneInteraction', 'Customized'), edgeM=NULL, layers=1, 
								filterGraphIDs=NULL, filterLayer=0, UP=TRUE, directed=FALSE, verbose=TRUE) {
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
			if (idType %in% c('GeneInteraction', 'Customized')) temp <- getSingleLayerGraphIDs(tempGraphIDs, edgeM, filterGraphIDs=tempFilter, directed=directed, UP=UP)
			else stop('The given idType is not supported!')
		}
		
		if (!directed && !is.null(temp)) {
			temp <- .removeDuplicatedEdge(temp)
			if (length(layersLength) > 0)  {
				l <- layersLength[length(layersLength)]
				tempL <- multiLayerGraphIDs[(length(multiLayerGraphIDs)-l+1):length(multiLayerGraphIDs)]
				tempLM <- .list2matrix(tempL, verbose=FALSE)
				tempCM <- .list2matrix(temp, verbose=FALSE)
				tempDLM <- rbind(tempLM, cbind(tempCM[,2], tempCM[,1]))
				tempCM <- tempCM[!duplicated(tempDLM)[-(1:dim(tempLM)[1])], ]
				if (!is.matrix(tempCM)) {
					tempCM <- matrix(tempCM, ncol=2)
				}
  				temp <- .matrix2list(tempCM, verbose=FALSE)
			}
			
		}
		
		
		if (!is.null(filterGraphIDs) && !is.null(temp)) {
			tempCM <- .list2matrix(temp, verbose=FALSE)
			tempSCNodes <- as.character(tempCM[tempCM[,1] == tempCM[,2], 1])
			tempNCNodes <- as.character(tempCM[!(tempCM[,1] == tempCM[,2]), ])
			tempSCNodes <- tempSCNodes[!(tempSCNodes %in% c(filterGraphIDs, tempNCNodes))]
			if (length(tempSCNodes) > 0) {
				temp <- temp[!(names(temp) %in% tempSCNodes)]
			}
		}
		
		if ((i == filterLayer) && (length(tempGraphIDs) > length(temp))) {
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
		
  		#if (is.null(multiLayerGraphIDs) | directed) multiLayerGraphIDs <- temp
		#else {
		#	tempMultiLayerGraphIDsM <- .list2matrix(multiLayerGraphIDs, verbose=FALSE)
		#	tempListM <- .list2matrix(temp, verbose=FALSE)
		#	tempIndex <- which(!duplicated(rbind(tempMultiLayerGraphIDsM, cbind(tempListM[,2], tempListM[,1])))) - nrow(tempMultiLayerGraphIDsM)
		#	multiLayerGraphIDs <- c(multiLayerGraphIDs, .matrix2list(tempListM[tempIndex[tempIndex > 0],], verbose=FALSE))
		#}
		
		tempGraphIDs <- unique(unlist(temp))
		tempGraphIDs <- tempGraphIDs[!(tempGraphIDs %in% unique(c(names(multiLayerGraphIDs), unlist(multiLayerGraphIDs), names(temp))))]
		multiLayerGraphIDs <- c(multiLayerGraphIDs, temp)
		layersLength <- c(layersLength, length(temp))
		#parentsGenes <- names(multiLayerGraphIDs)
		
		#tempGraphIDs <- tempGraphIDs[!(tempGraphIDs %in% parentsGenes)]
		
		if (length(tempGraphIDs) == 0) break
		
	}
	if (i < layers) selfClosure <- TRUE
	else selfClosure <- FALSE
	
	
	
#	if (!directed) {
#		removeDuplicates <- function(x,y) {
#			searchName <- unique(names(x))
#			searchIndex <- which(x < searchName)
#			otherIndex <- which(x >= searchName) 
#			if (length(searchIndex) != 0) {
#			    temp <- x[searchIndex]
#				#others <- x[which(x >= searchName)]
#				checkDuplicates <- function(p, q, r) {
#					if (r %in% q[[p]]) return(FALSE)
#					else return(TRUE)
#				}
#				return(x[sort(c(searchIndex[sapply(temp, checkDuplicates, y, searchName)], otherIndex))])
#			} else return(x)
#		}
#		
#		tempMatrix <- .list2matrix(multiLayerGraphIDs, verbose=FALSE)
#		rownames(tempMatrix) <- tempMatrix[,1]
#		tempList <- .matrix2list(tempMatrix, verbose=FALSE)
#		tempList <- lapply(tempList, removeDuplicates, tempList)
#		multiLayerGraphIDs <- tempList[which(sapply(tempList, length) > 0)]
#	}
	
	return(c(list(selfClosure), list(layersLength=layersLength), multiLayerGraphIDs))
}
