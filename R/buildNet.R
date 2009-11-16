`buildNet` <- 
function(graphIDs, idType=c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'GeneInteraction', 'Customized'), edgeM=NULL, layers=1, filterGraphIDs=NULL, filterLayer=0,
 					annLib=c('org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db', 'org.Dm.eg.db', 'customized'), output=c('interactive', 'fixed'), netMode=c('layer', 'connection'),
					vertexSize = NULL, edgeColor = NULL, colorMap=NULL, zeroColorIndex=NULL, matchMode=c('absolute', 'relative'),  label=TRUE, steric=FALSE, 
					directed=FALSE, direction=c('up', 'down', 'both'), showModeForNodes=c('nodes', 'filters'), verbose=TRUE, readable=TRUE, labelSize=1, labelColor='#666666',  ...) {
	netMode <- match.arg(netMode)
	annLib <- match.arg(annLib)
	idType <- match.arg(idType)
	output <- match.arg(output)
	matchMode <- match.arg(matchMode)
	showModeForNodes <- match.arg(showModeForNodes)
	direction <- match.arg(direction)
	if (is.null(filterGraphIDs)) filterLayer <- 0
	else {
		if (length(filterLayer) > 1) stop('filterLayer should be one integer!')
		else {
			if (filterLayer > layers) stop('filterLayer should not be more than total layers!')
		}
	}

	if (is.null(filterGraphIDs)) {
		filter <- FALSE
		if (filterLayer > 0) stop('Filter IDs are not available!')
	}
	else {
		filter <- TRUE
		if (!(is.vector(filterGraphIDs)) & !(is.data.frame(filterGraphIDs)) & !(is.matrix(filterGraphIDs))) stop('Filter IDs can not be recognized!')
		else {
			filterGraphIDs <- as.matrix(filterGraphIDs)
			if (dim(filterGraphIDs)[2] > 3) stop('Each filter ID has more than two specified value!')
			else {
				filterGraphIDs <- as.matrix(filterGraphIDs[!duplicated(filterGraphIDs),])
                if (length(filterGraphIDs[,1]) > length(unique(filterGraphIDs[,1]))) stop('One or more filter IDs have more than one specified values!')
				fGraphIDs <- matrix(rep('#000000', length(filterGraphIDs[,1])), ncol=1)
				rownames(fGraphIDs) <- filterGraphIDs[,1]
			}
		}
	}
	graphIDs <- unique(as.character(unlist(graphIDs)))
	if (NA %in% graphIDs) {
		print('Warning: IDs contain NA! Removing NA from IDs')
		root <- graphIDs[!(graphIDs %in% NA)]
	}
	if (netMode == 'layer') {
		if (is.null(vertexSize)) vsize <- 5 * 1.5^(1.2 * c(layers:0))
		else {
			if (length(vertexSize) != (layers+1)) stop('vertexSize contains different layers!')
			if (is.numeric(vertexSize) & is.vector(vertexSize)) vsize <- vertexSize
			else stop('vertexSize is not valid!')
		}
		if (is.null(edgeColor)) ecolor <- c('#0077ff', '#ff9900', '#660033', '#669900', '#cc00ff', '#339966', '#9999cc', '#ccff33', '#6600ff', '#ff3399', rep('#ff3399', 20)) 
		else {
			if (length(edgeColor) != layers) stop('edgeColor contains different layers!')
			if (is.character(edgeColor) & is.vector(edgeColor)) ecolor <- edgeColor
			else stop('edgeColor is not valid!')
		}
	} else {
		vsize <- 5 * 2^(1.1 * c(2,rep(1, layers)))
		ecolor <- c('#6600ff', rep('#ff9900', (layers-1)))
	}
	
	extendIA <- function(IA) {
		if (dim(IA)[2] > 2) tempM <- cbind(IA[,2], IA[,1], IA[,3:dim(IA)[2]])
		else tempM <- cbind(IA[,2], IA[,1])
		colnames(tempM) <- colnames(IA)
		tempM <- rbind(IA, tempM)[,1:2]
		return(tempM[!(duplicated(tempM)),])
	}
	if (!all(idType %in% c('GO', 'GO.BP', 'GO.CC', 'GO.MF'))) {
		if (idType == 'GeneInteraction') {
			if (annLib == 'customized') {
				if (is.null(edgeM)) stop("Customized database is not available!")
			} else {
				# if DOLite provides ontology structure, here should be a flow control for edgeM.
				#load('~/Documents/projects/GeneInteraction/HsIALite.rda')
				#load('~/Documents/projects/GeneInteraction/MmIALite.rda')
				#load('~/Documents/projects/GeneInteraction/RnIALite.rda')
				#load('~/Documents/projects/GeneInteraction/DmIALite.rda')
				switch(annLib,
					'org.Hs.eg.db'=data('HsIALite', package='GeneAnswers'),
					'org.Mm.eg.db'=data('MmIALite', package='GeneAnswers'),
					'org.Rn.eg.db'=data('RnIALite', package='GeneAnswers'),
					'org.Dm.eg.db'=data('DmIALite', package='GeneAnswers'))
				edgeM <- switch(annLib,
					'org.Hs.eg.db'=extendIA(HsIALite),
			   		'org.Mm.eg.db'=extendIA(MmIALite),
					'org.Rn.eg.db'=extendIA(RnIALite),
					'org.Dm.eg.db'=extendIA(DmIALite))
			}
		} else {
			if (is.null(edgeM)) stop("Customized database is not available!")
		}
	}
	
	if (verbose) print('Building graph structure ...')
	
	inputList <- c()
	UP = FALSE
	if (filter) {
		filterIDs=filterGraphIDs[,1]
		filterLayer=filterLayer
	} else {
		filterIDs=NULL
		filterLayer=0
	}
	
	if (directed) {
		if (direction %in% c('down', 'both')) {
			inputList <- getMultiLayerGraphIDs(graphIDs, idType=idType, edgeM=edgeM, layers=layers, filterGraphIDs=filterIDs, filterLayer=filterLayer, UP=FALSE, directed=directed, verbose=verbose)
		}
		if (direction %in% c('up', 'both')) {
			upInputList <- getMultiLayerGraphIDs(graphIDs, idType=idType, edgeM=edgeM, layers=layers, filterGraphIDs=filterIDs, filterLayer=filterLayer, UP=TRUE, directed=directed, verbose=verbose)
			if (is.null(inputList[[2]])) {
				inputList <- upInputList
				UP=TRUE
			} else {
				if (!is.null(upInputList[[2]])) {
					downInputList <- inputList
					inputList[[1]] <- TRUE
					inputList[[2]] <- c(upInputList[[2]], downInputList[[2]])
					tempUpM <- .list2matrix(upInputList[-1:-2], verbose=FALSE)
					inputList <- c(inputList[1:2], .matrix2list(cbind(tempUpM[,2], tempUpM[,1]),  verbose=FALSE), downInputList[-1:-2])
				}
			}
		}
	} else {
		inputList <- getMultiLayerGraphIDs(graphIDs, idType=idType, edgeM=edgeM, layers=layers, filterGraphIDs=filterGraphIDs, filterLayer=filterLayer, UP=TRUE, directed=directed, verbose=verbose)
	}
	

	selfClosure <- inputList[[1]]
	inputList <- inputList[2:length(inputList)] 
	returnResult <- inputList[2:length(inputList)]
    
	if (is.null(unlist(inputList[2:length(inputList)]))) {
		print('All nodes are not connected! Aborting ...')
		return(NULL)
	}
	
	temp <- inputList[[1]]
	if (readable) {
		if (idType %in% c('GeneInteraction')) {
			if (annLib != 'customized') {
				require(annLib, character.only=TRUE)
		    	inputList <- c(inputList[1], lapply(inputList[2:length(inputList)], getSymbols, annLib))
				names(inputList)[2:length(inputList)] <- getSymbols(names(inputList)[2:length(inputList)],annLib)
				if (filter) {
					rownames(fGraphIDs) <- getSymbols(rownames(fGraphIDs), annLib)
					filterGraphIDs[,1] <-  getSymbols(filterGraphIDs[,1], annLib)
				}
				graphIDs <- getSymbols(graphIDs, annLib)
			}
		} else {
			inputList <- c(inputList[1], lapply(inputList[2:length(inputList)], getCategoryTerms, idType, ...))
			names(inputList)[2:length(inputList)] <- getCategoryTerms(names(inputList)[2:length(inputList)], idType, ...)
			if (filter) {
				rownames(fGraphIDs) <- getCategoryTerms(rownames(fGraphIDs), idType, ...)
				filterGraphIDs[,1] <-  getCategoryTerms(filterGraphIDs[,1], idType, ...)
			}
			graphIDs <- getCategoryTerms(graphIDs, idType, ...)
		}
	} 
	
	
	tempG <- .list2graph(inputList[2:length(inputList)], directed=directed, verbose=verbose, reverse=UP)
	g <- tempG[[1]]
	if (verbose) {
		if (directed) {
			tempIn <- igraph::degree(g, mode='in')
			tempOut <- igraph::degree(g, mode='out')
			names(tempOut) <- c(1:length(tempOut))
			possibleRoots <- tempOut[which(tempIn == min(tempIn))]
	   		print(paste('For the given directed graph, the node ', paste(as.character(as.numeric(names(which(possibleRoots == max(possibleRoots)))) - 1), collapse=' or '),' might be the root.', sep=''))
	   	}
	   	else {
	   		print(paste('For the given undirected graph, the node ', (which.max(igraph::degree(g)) - 1),' might be the root.', sep=''))
	   	}
	}
	
	
	if (filter){
		if (showModeForNodes == 'nodes') {
			numberOfCol = dim(filterGraphIDs)[2]
			filterGraphIDs <- matrix(filterGraphIDs[filterGraphIDs[,1] %in% unique(c(names(tempG[[2]]), graphIDs)),], ncol=numberOfCol)
	        tempFGraphIDs <- fGraphIDs[rownames(fGraphIDs) %in% unique(c(names(tempG[[2]]), graphIDs)),]
			fGraphIDs <- matrix(tempFGraphIDs, ncol=1)
			rownames(fGraphIDs) <- names(tempFGraphIDs)
		}
		if (dim(filterGraphIDs)[2] >= 2) {
			if (all(!is.na(as.numeric(filterGraphIDs[,2])))) {
				if (length(unique(as.numeric(filterGraphIDs[,2]))) > 1) {
					if (is.null(colorMap))  {
						colorLevel <- 256
						#zeroColorIndex <- ceiling(colorLevel/2)
						if ((min(as.numeric(filterGraphIDs[,2])) > 0) | (max(as.numeric(filterGraphIDs[,2])) < 0)) {
							if (min(as.numeric(filterGraphIDs[,2])) > 0) {
								conceptCol <- colorRampPalette(c('#ffb7b7','#ff0000'))
							} else {
								conceptCol <- colorRampPalette(c('#00ff00','#b7ffb7'))
							}
							matchMode <- 'relative'
							colorMap <- conceptCol(colorLevel)
						} else {
							zeroColorIndex <- 1 + ceiling(abs(min(as.numeric(filterGraphIDs[,2]))) * (colorLevel-1) / (max(as.numeric(filterGraphIDs[,2]))-min(as.numeric(filterGraphIDs[,2]))))
							conceptCol <- colorRampPalette(c('#00ff00','#ffffff'))                  
							colorMap <- conceptCol(zeroColorIndex)
							conceptCol <- colorRampPalette(c('#ffffff','#ff0000'))
							colorMap <- c(colorMap, conceptCol(colorLevel-zeroColorIndex + 1))
						}
					}
					fGraphIDs[,1] <- .colorMatch(as.numeric(filterGraphIDs[,2]), colorMap=colorMap, matchMode=matchMode, zeroColorIndex=zeroColorIndex)
				} else {
					if (is.null(colorMap)) fGraphIDs[,1] <- '#000000'
					else fGraphIDs[,1] <- colorMap[1]
				}
			} else {
				stop('The 2nd column can not be recognized!')
			}
			if (dim(filterGraphIDs)[2] == 3) {
				if (all(!is.na(as.numeric(filterGraphIDs[,3])))) fGraphIDs <- cbind(fGraphIDs, as.numeric(filterGraphIDs[,3]))
				else stop('The 3nd column can not be recognized!')
			}
		}
	}

	if (netMode == 'layer') E(g)$color <- '#aaaaaa'
	else E(g)$color <- '#ff9900'
	E(g)$width <- 3
	V(g)$color <- "#ffffff"
	
	exclude <- function(x, y) {
		result <- list()
		if (length(x) == length(y)) {
			for (i in 1: length(x)) {
				result <- c(result, x[[i]][!(x[[i]] %in% y[[i]])])
			}
			return(result) 
		}
	} 
	if (selfClosure) {
		layerStr <- list(tempG[[2]][names(tempG[[2]]) %in% graphIDs])
		parentLayer <- igraph::neighborhood(g,0, nodes=tempG[[2]][names(tempG[[2]]) %in% graphIDs])
	}
	else {
		layerStr <- list(c(0:(temp[1]-1)))
		parentLayer <- igraph::neighborhood(g,0, nodes=0:(temp[1]-1))
	}
	for (i in 1:layers) {
		if (selfClosure) totalLayer <- igraph::neighborhood(g,i, nodes=tempG[[2]][names(tempG[[2]]) %in% graphIDs])
		else totalLayer <- igraph::neighborhood(g,i,nodes=0:(temp[1]-1))
		currentLayer <- exclude(totalLayer, parentLayer)
		currentLayer <- lapply(currentLayer, function(x, y) {return(x[!(x %in% y)])}, unique(unlist(layerStr)))
		tempStr <- unique(unlist(currentLayer))
		if (all(tempStr %in% layerStr)) {
			i <- i - 1;
			break
		}
		layerStr <- c(layerStr, list(tempStr))
		parentLayer <- totalLayer
	}
	
    layers <- min(layers, i)
	V(g)[tempG[[2]][names(tempG[[2]]) %in% graphIDs]]$size <- vsize[1]
	if (layers > 0) {
		for (i in 1:layers) {
			E(g)[V(g)[layerStr[[i]]] %--% V(g)[layerStr[[i+1]]]]$color <- ecolor[i]
			V(g)[layerStr[[i+1]]]$size <- vsize[i+1]
		}
	}
	if (filter) V(g)[tempG[[2]][names(tempG[[2]]) %in% rownames(fGraphIDs)]]$color <- fGraphIDs[names(tempG[[2]][names(tempG[[2]]) %in% rownames(fGraphIDs)]),1]

	if (label) verticesLabel <- names(tempG[[2]])
	else {
		verticesLabel <- rep('', length(tempG[[2]]))
		verticesLabel[(tempG[[2]][names(tempG[[2]]) %in% graphIDs] + 1)] <- names(tempG[[2]])[(tempG[[2]][names(tempG[[2]]) %in% graphIDs] + 1)]
	}	
	
	E(g)[tempG[[2]][names(tempG[[2]]) %in% graphIDs] %--% tempG[[2]][names(tempG[[2]]) %in% graphIDs]]$color <- ecolor[1]
	
	if (!all(graphIDs %in% names(tempG[[2]]))) {
		singleNodes <- graphIDs[!(graphIDs %in% names(tempG[[2]]))]
		g <- add.vertices(g, length(singleNodes))
		tempNodes <- length(tempG[[2]]) + c(1:length(singleNodes)) - 1
		names(tempNodes) <- singleNodes
		tempG[[2]] <- c(tempG[[2]], tempNodes)
		V(g)[tempNodes]$size <- vsize[1]
		V(g)[tempNodes]$color <- '#ffffff'
		if (filter) V(g)[tempNodes[rownames(fGraphIDs)[rownames(fGraphIDs) %in% singleNodes]]]$color <- fGraphIDs[rownames(fGraphIDs) %in% singleNodes,1]
		verticesLabel <- c(verticesLabel, singleNodes)
	}
	
	if (filter) {
		if (dim(fGraphIDs)[2] == 2) {
			fNodes <- intersect(names(tempG[[2]]), rownames(fGraphIDs))
			V(g)[tempG[[2]][fNodes]]$size <- 5*sqrt(abs(as.numeric(fGraphIDs[fNodes, 2])))
		}
	}
	
	nodeFrameColor <- rep('black', length(V(g)))
	nodeFrameColor[tempG[[2]][names(tempG[[2]]) %in% graphIDs]+1] <- 'yellow'
	V(g)$label = verticesLabel
	V(g)$label.color = labelColor
	V(g)$frame.color = nodeFrameColor 
	if (output == 'interactive') {
		tkplot(g, vertex.label.font=2, vertex.label=verticesLabel, vertex.label.color=labelColor, vertex.label.cex=labelSize, vertex.frame.color=nodeFrameColor, edge.loop.angle=-90, layout=layout.fruchterman.reingold.grid, canvas.width=800, canvas.height=600)
	} else {
		plot(g, vertex.label.font=2, vertex.label=verticesLabel, vertex.label.color=labelColor, vertex.label.cex=labelSize, vertex.frame.color=nodeFrameColor, edge.loop.angle=-90, layout=layout.fruchterman.reingold.grid)
	}
		
	if (steric) {
		require(rgl)
		bg3d('#555555')
		coordinateM <- layout.fruchterman.reingold(g, dim=3)
		rglplot(g, vertex.label.font=2, vertex.size=(V(g)$size)/2, vertex.color=V(g)$color, vertex.label=verticesLabel, vertex.label.color='#aaaaaa', vertex.label.cex=4, vertex.label.dist=0.35, edge.width=1, layout=coordinateM)
	}
	return(invisible(c(list(g), returnResult)))
}
