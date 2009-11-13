.onLoad <- function(libname, pkgname) {
	if(.Platform$OS.type == "windows" && require(Biobase) && interactive()
        && .Platform$GUI ==  "Rgui") {
        addVigs2WinMenu("GeneAnswers")
    }

    #if (.Platform$OS.type == "unix" && (.Platform$GUI %in% c("X11", "Tk", "GNOME"))) {
	#	quartz <- function(...) X11(...)
	#}
	#
	#if (.Platform$OS.type == "windows" && (.Platform$GUI %in% c("Rgui", "Rterm"))) {
	#	quartz <- function(...) windows(...)
	#}
	 
	data('DOLite', package='GeneAnswers')
	data('DOLiteTerm', package='GeneAnswers')
	data('HsIALite', package='GeneAnswers')
	data('MmIALite', package='GeneAnswers')
	data('RnIALite', package='GeneAnswers')
	data('DmIALite', package='GeneAnswers')
	
}

.filterGOIDs <- function(GOCategory=c('ALL', 'BP', 'CC', 'MF'), level=level) {
	require(GO.db)
	GOCategory <- match.arg(GOCategory)
	filterGOIDs <- c()
	if (level > 0) {
		if ((GOCategory == 'ALL') | (GOCategory == 'BP')) filterGOIDs <- unique(c(filterGOIDs, 'GO:0008150'))
		if ((GOCategory == 'ALL') | (GOCategory == 'CC')) filterGOIDs <- unique(c(filterGOIDs, 'GO:0005575'))
		if ((GOCategory == 'ALL') | (GOCategory == 'MF')) filterGOIDs <- unique(c(filterGOIDs, 'GO:0003674'))
		if (level > 1) {
			for (i in 2:level) {                                                                              
				if ((GOCategory == 'ALL') | (GOCategory == 'BP')) filterGOIDs <- unique(c(filterGOIDs, unlist(lookUp(filterGOIDs, "GO", "BPCHILDREN"))))
				if ((GOCategory == 'ALL') | (GOCategory == 'CC')) filterGOIDs <- unique(c(filterGOIDs, unlist(lookUp(filterGOIDs, "GO", "CCCHILDREN"))))
				if ((GOCategory == 'ALL') | (GOCategory == 'MF')) filterGOIDs <- unique(c(filterGOIDs, unlist(lookUp(filterGOIDs, "GO", "MFCHILDREN"))))
				filterGOIDs <- unique(filterGOIDs[!(filterGOIDs %in% NA)]) 
			}
		}
	}
	return(filterGOIDs)
}

#.makeGraph <- function (edgesMatrix=NULL, openFile=FALSE, fileName=NULL, saveFile = FALSE) {
#	if (openFile) edgesMatrix <- as.matrix(read.table(fileName, sep='\t'))
#	originalNodes <- unique(as.list(edgesMatrix))
#	nodes <- rep(0:(length(originalNodes)-1))
#	names(nodes) <- originalNodes
#	newEdges <- matrix(nodes[as.character(edgesMatrix[,1])])
#	newEdges <- cbind(newEdges, nodes[as.character(edgesMatrix[,2])]) 
#	if (saveFile) {
#		write.table(newEdges, file='graph.txt', sep=' ', row.names=F, col.names=F)
#		return(nodes)
#	} else {
#		return(c(list(newEdges), list(nodes)))
#	}
#}
#
#.list2graph <- function (inputList, ...) {
#	column1 <- c()
#	if (length(inputList) >0) {
#		for (i in 1:length(inputList)) {
#			column1 <- c(column1, rep(-i, length(inputList[[i]])))
#		}
#		edgeMatrix <- matrix(column1)
#		edgeMatrix <- cbind(edgeMatrix, unlist(inputList))
#		return(.makeGraph(edgesMatrix=edgeMatrix, ...))
#	} else {
#		stop('list is empty!')
#	}
#}

.makeGraph <- function (edgesMatrix=NULL, directed=FALSE, openFile=FALSE, fileName=NULL, reverse=FALSE) {
	require(igraph)
	if (openFile) {
		edgesMatrix <- as.matrix(read.table(fileName, sep='\t'))
		if (dim(edgesMatrix)[2] != 2) stop('Edge matrix should contain only 2 columns!')
	}
	originalNodes <- unique(as.vector(edgesMatrix))
	nodes <- rep(0:(length(originalNodes)-1))
	names(nodes) <- as.character(originalNodes)
	newEdges <- matrix(nodes[as.character(edgesMatrix[,1])])
	newEdges <- cbind(newEdges, nodes[as.character(edgesMatrix[,2])])
	if (reverse) {
		newEdges <- cbind(newEdges[,2], newEdges[,1])
	}
	g <- graph.edgelist(newEdges, directed=directed) 
	return(c(list(g), list(nodes)))
}


.list2graph <- function (inputList, verbose=TRUE, directed=FALSE, ...) {
	if (!(all(!(is.null(names(inputList)))))) stop('Names of inputlist shoule not contain NULL!')
	if (verbose) print('Removing element containing NULL')
	inputList <- inputList[which(sapply(inputList, length) != 0)]
	column1 <- c()
	if (verbose) print('Removing NA')
	removeNA <- function(x) { x <- x[!(x %in% NA)]; return(as.character(unique(x)))}
	tempList <- lapply(inputList, removeNA) 
	tempList <- tempList[which(sapply(tempList, length) != 0)]
	if (length(tempList) > 0) {
		for (i in 1:length(tempList)) {
			column1 <- c(column1, rep(names(tempList)[i], length(tempList[[i]])))
		}
		edgeMatrix <- matrix(column1)
		edgeMatrix <- cbind(edgeMatrix, unlist(tempList))
		edgeMatrix <- edgeMatrix[!(duplicated(edgeMatrix)),]
		if (!is.matrix(edgeMatrix)) edgeMatrix <- matrix(edgeMatrix, ncol=2)
		return(.makeGraph(edgesMatrix=edgeMatrix, directed=directed, ...))
	} else {
		stop('list is empty!')
	}
}

.drawTable <- function(dataMatrix, mar=c(1,1,5,8), addRowLabel=TRUE, cex.axis=c(1.1, 0.9), ...) {
	if (is.null(rownames(dataMatrix)) | is.null(colnames(dataMatrix))) stop('rownames and/or colnames of input matrix are missing!')
	oldsetting <- par('mar'=mar) 
	#image(x=1:ncol(dataMatrix), y=1:nrow(dataMatrix), -t(dataMatrix), ylim=c((nrow(dataMatrix)+1),0), col=c('white','white'), axes=F, xlab='', ylab='', )
	image(x=1:ncol(dataMatrix), y=1:nrow(dataMatrix), -t(dataMatrix), col=c('white','white'), axes=F, xlab='', ylab='', )
	if(nrow(dataMatrix) < 120) abline(h=c(0:nrow(dataMatrix))+0.5, col='#bbbbbb', lty = "solid", lwd=2)
	if (addRowLabel) axis(2,at=1:nrow(dataMatrix), labels=rownames(dataMatrix), tick=F,las=2, cex.axis=cex.axis[length(cex.axis)])
	axis(3,at=1:ncol(dataMatrix), labels=colnames(dataMatrix), tick=F,las=2, cex.axis=cex.axis[1])
	OnesPos <- which(dataMatrix== 1, arr.ind=TRUE)
	points(OnesPos[,2], OnesPos[,1], pch=19, col='black', cex=18/nrow(dataMatrix))
	par(oldsetting)
}

# dataMatrix: a 2-dimensional numeric matrix
# sortBy:	determine whether to sort the dataMatrix by row, column, both row and column or none of them,
#			The sorting is calculated by function isoMDS
# rotate: whether to rotate the data matrix or not
# standardize: whether to standardize the matrix by row
# colorMap: color map of the heatmap
# mar: margin of the plot
# lib: annotation library for converting probe Id (row) to gene symbol. If it is NULL, it will do nothing
# cex.axis: the character size of row and column labels
# maxQuantile: Set the satuation values based on quantile (for the purpose of better visuallizaiton by removing possible outliers)
# maxVal: Set the satuation values based on absolute values (for the purpose of better visuallizaiton by removing possible outliers) 
#		Both maxQuantile and maxVal can be a vector of lenght 2, which is the down and up limit of the range
# addRowLabel: whether to add row lables or not
# rm.unannotate: remove un-annotated probes when lib is provided.
# reverseSort: reverse sorting of the rows or columns
# colorBar: determine whether plot the color bar in the separated window (In the current implementation, users have to the adjust the window to get a better view)
# colorBarLabel: labels of color bar
.heatmap.mds <- function(dataMatrix, sortBy=c('row', 'column', 'both', 'none'), rotate=FALSE, standardize=TRUE, colorMap='RdYlBu', mar=c(1,1,5,8),  
	cex.axis=c(0.9, 0.9), maxQuantile=0.99, maxVal=NULL, symmetry=FALSE, addRowLabel=TRUE, rm.unannotate=TRUE, reverseSort=c('none', 'row', 'column', 'both'), 
	colorBar=FALSE, colorBarLabel=NULL, mapType=c('table', 'heatmap'),  ...)
{
	if (colorMap[1] == 'GBR') {
		library(Heatplus)
		colorMap <- rev(RGBColVec(256))
	} else {
		require(RColorBrewer)
	}
	require(MASS)
	sortBy <- match.arg(sortBy)
	reverseSort <- match.arg(reverseSort)
	sort <- TRUE
	if (sortBy == 'none') sort <- FALSE
	dataRange <- range(dataMatrix)
	if (symmetry) {
		temp <- max(abs(dataRange))
		dataRange <- c(-temp, temp)
	} 
	
	if (!is.null(maxQuantile)) {
		if (length(maxQuantile) == 1) maxQuantile <- c(1-maxQuantile, maxQuantile)
		maxQuantile <- sort(maxQuantile)
		if (any(maxQuantile > 1)) maxQuantile[maxQuantile > 1] <- 1
		if (any(maxQuantile < 0)) maxQuantile[maxQuantile < 0] <- 0
		maxQVal <- quantile(dataMatrix, maxQuantile)
		if (!is.null(maxVal)) {
			maxVal <- sort(maxVal)
			if (maxVal[1] < maxQVal[1]) maxVal[1] <- maxQVal[1]
			if (maxVal[2] > maxQVal[2]) maxVal[2] <- maxQVal[2]
		} else {
			maxVal <- maxQVal
		}
	}

	if (!is.null(maxVal)) {
		if (length(maxVal) == 1) maxVal <- c(-maxVal, maxVal)
		maxVal <- sort(maxVal)
		maxVal <- round(maxVal, 2)
		dataMatrix[dataMatrix < maxVal[1]] <- maxVal[1]
		dataMatrix[dataMatrix > maxVal[2]] <- maxVal[2]
	}	
	
	ind.order <- list(row=1:nrow(dataMatrix), column=1:ncol(dataMatrix))
	## get the gene symbol
	rowLabels <- rownames(dataMatrix)
	if (length(rowLabels) < 2) {
		warning('The dataMatix should have more than 2 rows!')
		return(NULL)
	}

	if (rotate) {
		colLabels <- rowLabels
		rowLabels <- colnames(dataMatrix)
	} else {
		colLabels <- colnames(dataMatrix)
	}

	## Standardize the gene expression
	if (standardize) {
		relative <- TRUE
		dataMatrix <- t(scale(t(dataMatrix)))
		if (any(is.na(dataMatrix))) {
			mm <- mean(dataMatrix, na.rm=T)
			dataMatrix[is.na(dataMatrix)] <- mm
		}
	} else {
		relative <- FALSE
	}

	seqColor <- c('BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral')
	if (length(colorMap) == 1) {
		if (colorMap %in% seqColor) colorMap <- colorRampPalette(rev(brewer.pal(11, colorMap)))(256)
	}

	## reorder the genes based on the MDS distance
	# start with a distance matrix
	if (sort) {
		if (sortBy == 'row' || sortBy == 'both') {
			distance.matrix <- as.matrix(dist(dataMatrix, method='euclidean', upper=T))
			distance.matrix[distance.matrix < 1e-6] <- 1e-6
			diag(distance.matrix) <- 0
			result.mds <- isoMDS(distance.matrix, k=1)
			# result.mds <- sammon(distance.matrix)
			i.ordered <- order(result.mds$points)
			if (reverseSort == 'row' || reverseSort == 'both') i.ordered <- rev(i.ordered)
			if (rotate) {
				dataMatrix <- t(dataMatrix[i.ordered, ])
				colLabels <- colLabels[i.ordered]
			} else {
				dataMatrix <- dataMatrix[i.ordered, ]
				rowLabels <- rowLabels[i.ordered]			
			}
			ind.order$row <- i.ordered
		} 
		if (sortBy == 'column' || sortBy == 'both') {
			distance.matrix <- as.matrix(dist(t(dataMatrix), method='euclidean', upper=T))
			distance.matrix[distance.matrix < 1e-6] <- 1e-6
			diag(distance.matrix) <- 0
			result.mds <- isoMDS(distance.matrix, k=1 )
			# result.mds <- sammon(distance.matrix)
			i.ordered <- order(result.mds$points)
			if (reverseSort == 'column' || reverseSort == 'both') i.ordered <- rev(i.ordered)
			if (rotate) {
				dataMatrix <- t(dataMatrix[, i.ordered])
				rowLabels <- rowLabels[i.ordered]
			} else {
				dataMatrix <- dataMatrix[, i.ordered]
				colLabels <- colLabels[i.ordered]
			}
			ind.order$column <- i.ordered
		}
	} else {
		if (rotate) {
			dataMatrix <- t(dataMatrix)
		}
	}
	
	
	if (mapType == 'table') .drawTable(dataMatrix, addRowLabel=addRowLabel, mar=mar, cex.axis=cex.axis, ...)
	
	else {
		oldsetting <- par('mar'=mar)
		if (symmetry) {
			maxV <- max(abs(range(dataMatrix)))
			maxVal <- c(-maxV, maxV)
			image(1:ncol(dataMatrix), 1:nrow(dataMatrix), t(dataMatrix), col=colorMap, axes=F, xlab='', ylab='', zlim=maxVal, ...)
		}
		else image(1:ncol(dataMatrix), 1:nrow(dataMatrix), t(dataMatrix), col=colorMap, axes=F, xlab='', ylab='', ...)  
		#image(x=1:ncol(inputM), y=1:nrow(inputM), -t(inputM), ylim=c((nrow(inputM)+1),0), col=c('white','white'), axes=F, xlab='', ylab='', )
		axis(3, at=1:ncol(dataMatrix), labels=colLabels, tick=F, las=2, cex.axis=cex.axis[1])
		if (addRowLabel) {
			axis(4, at=1:nrow(dataMatrix), labels=rowLabels, tick=F, las=2, cex.axis=cex.axis[length(cex.axis)])	
		} else {
			axis(4, at=1:nrow(dataMatrix), labels=rep('', length(rowLabels)), tick=F, las=2, cex.axis=cex.axis[length(cex.axis)])
		}
		box()
		
		if (colorBar) {
			deviceNo <- dev.cur()	
			x11()
			par(mar=c(2,2,3,5))
			color <- seq(maxVal[1], maxVal[2], length=length(colorMap))
			colorBar <- t(matrix(rep(color,2), ncol=2))
			image(1:nrow(colorBar), 1:ncol(colorBar), colorBar, col=colorMap, axes=F, xlab="")
			if (relative) {
				if (is.null(colorBarLabel)) {
					colorBarLabel <- c('low', 'high')
				}
				axis(4, at = c(1.5, length(colorMap)-0.5), labels=colorBarLabel, tick=F)			
			} else {
				if (is.null(colorBarLabel)) {
					colorBarLabel <- seq(dataRange[1], dataRange[2], length=7)
				}
				if (is.numeric(colorBarLabel)) 
					colorBarLabel <- round(colorBarLabel,2)
				if (!is.null(maxVal)) {
					if (colorBarLabel[1] == -round(maxVal[1], 2)) {
						colorBarLabel[1] <- paste('-', round(maxVal, 2), '<')
					}
					if (colorBarLabel[length(colorBarLabel)] == round(maxVal[1], 2)) {
						colorBarLabel[length(colorBarLabel)] <- paste('>', round(maxVal, 2))
					}
				}
				axis(4, at = seq(1, length(colorMap), length=length(colorBarLabel)), labels=colorBarLabel)		
			}
			box()
			dev.set(deviceNo)
		}

		par(oldsetting)
	}
	return(invisible(ind.order))
}

.hyperGTest <- function(inputVector, categoryList, ...) {
	if (is.vector(inputVector) & is.character(inputVector)) {
		testResult <- lapply(categoryList, .overReprsTest, inputVector, ...)
		hyperGTest <- as.data.frame(t(matrix(unlist(testResult), ncol=length(testResult), nrow=length(testResult[[1]]))))
		rownames(hyperGTest) <- names(testResult)
		colnames(hyperGTest) <- names(testResult[[1]])
		return(cbind(hyperGTest, 'fdr p value'=p.adjust(hyperGTest[,'p value'], method='fdr')))
	} else {
		stop('Input should be a character vector! ') 
	}
}

.overReprsTest <- function(categoryIDs, observedIDs, totalNGenes = 30000) {
	if ((length(observedIDs) > 0) & (length(categoryIDs) > 0)) {
		observedIDs <- unique(unlist(observedIDs))
		categoryIDs <- unique(unlist(categoryIDs))
		inCategoryIDs <- observedIDs[observedIDs %in% categoryIDs]
		perInObserved <- length(inCategoryIDs)/length(observedIDs)
		perInGenome <- length(categoryIDs)/totalNGenes
		foldOverRepresents <- perInObserved/perInGenome
		odds_ratio <- (length(inCategoryIDs) * (totalNGenes - length(observedIDs) - length(categoryIDs) + length(inCategoryIDs))) / 
						((length(categoryIDs) - length(inCategoryIDs)) * (length(observedIDs) - length(inCategoryIDs)))
		output <- c(length(inCategoryIDs), perInObserved, perInGenome, foldOverRepresents, odds_ratio, 
						phyper(length(inCategoryIDs) - 1L, length(categoryIDs), totalNGenes-length(categoryIDs), length(observedIDs), lower.tail = FALSE))
	} else {
		output <- c(0, NA, NA, NA, NA, 1)
	}
	names(output) <- c('genes in Category', 'percent in the observed List', 'percent in the genome', 'fold of overrepresents', 'odds ratio', 'p value')
	return(output)
}

.searchEntrezTag<- function(tag, species='Homo sapiens') {
  unlist(lapply(tag, .searchEntrezTerm, species=species))
}

.searchEntrezTerm <- function(term, baseUrl="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/", species=c('human', 'rat', 'mouse', 'fly')) {
	require(XML)
	species <- match.arg(species)
	species = switch(species,
		'human'='Homo sapiens',
  		'rat'='Rattus norvegicus',
  		'mouse'='Mus musculus',
  		'fly'='Drosophila melanogaster')
	checkHttpSyntax <- function(x, type=c('quote', 'plus')) {
		if (length(grep(pattern=" ", x)) != 0) {
		   if (type == 'quote') x <- paste('"', x, '"', sep='')
		   else x <- gsub(' ', '+', x)
		}
		return(x)
	}
    # QC: make sure the baseUrl is all right.
    if (is.null(baseUrl)) {
        stop("Need to define the URL of the Pubmed service!")
    } 
	
	if (length(term) > 0) {
		# replace space for Http syntax 
		term <- checkHttpSyntax(term, type='quote')
		species <- checkHttpSyntax(species, type='plus')
	  	# Get the query string ready. This string should be in the Pubmed 
	  	# syntax. The Pubmed syntax is documented at
	  	# http://eutils.ncbi.nlm.nih.gov/entrez/query/static/esearch_help.html
	    query <- URLencode(paste(baseUrl, "esearch.fcgi?", "db=gene&term=", term, "+AND+",
					   species, "[organism]&retmode=xml&retmax=500000", sep=""))
	} else {
		species <- checkHttpSyntax(species, type='plus')
	  	# Get the query string ready. This string should be in the Pubmed 
	  	# syntax. The Pubmed syntax is documented at
	  	# http://eutils.ncbi.nlm.nih.gov/entrez/query/static/esearch_help.html
	    query <- URLencode(paste(baseUrl, "esearch.fcgi?", "db=gene&term=", 
					   species, "[organism]&retmode=xml&retmax=500000", sep=""))
	} 
	print(paste('search link:', query))  
    # parse resulting XML into a tree to be returned to user
    result.xml <- try(xmlTreeParse(file=query, isURL=TRUE))
    #entrezIDs <- as.numeric(xmlSApply(xmlRoot(result.xml)[	["IdList"]], xmlValue))
	entrezIDs <- xmlSApply(xmlRoot(result.xml)[["IdList"]], xmlValue)
    names(entrezIDs) <- rep(term, length(entrezIDs))
    return(entrezIDs)
}

.getGOTerms <- function(GOIDs) {
	temp <- getGOTerm(GOIDs)
	names(temp) <- NULL
	return(unlist(temp)[GOIDs])
}
 
.colorMatch <- function(values, colorMap=NULL, matchMode=c('absolute', 'relative'), zeroColorIndex=NULL) {
	matchMode <- match.arg(matchMode)
	tempColor <- c()

	colorMapping <- function(inputValues, inputColors) {
		minValue <- min(inputValues)
		maxValue <- max(inputValues)
		outputColors <- inputColors[1 + floor((inputValues-minValue) * (length(inputColors))/((maxValue-minValue)*1.001))]
		names(outputColors) <- as.character(inputValues)
		return(outputColors)
	}

	if (matchMode == 'absolute') {
		if (is.numeric(zeroColorIndex) & (length(zeroColorIndex) == 1)) {
			if ((zeroColorIndex < length(colorMap)) & (zeroColorIndex > 1)) {
				tempColor <- c(tempColor, colorMapping(c(0, values[values >= 0]), colorMap[zeroColorIndex:length(colorMap)])[-1])
				tempColor <- c(tempColor, colorMapping(c(0, values[values < 0]), colorMap[1:(zeroColorIndex-1)])[-1])
				outColor <- tempColor[as.character(values)]
			} else stop('zeroColorIndex is out of colorMap index! Aborting ...')
		} else stop('zeroColorIndex is not ONE NUMBER! Aborting ...')
	} else {
		outColor <- colorMapping(values, colorMap)
	}
	return(outColor) 
}

.getSingleGraphIDs <- function(graphID, edgeM, filterGraphIDs=NULL, UP=TRUE) {
	if (!is.character(graphID)) graphID <- as.character(graphID)
	options(warn=-1)
	if (!is.na(as.character(graphID))) {
		if (UP) temp <- as.character(edgeM[as.character(edgeM[,1]) %in% graphID,2])
		else temp <- as.character(edgeM[as.character(edgeM[,2]) %in% graphID,1])
		if (!is.null(filterGraphIDs)) temp <- intersect(temp, as.character(filterGraphIDs))
		options(warn=0)
		if (length(temp) > 0) return(unique(temp[!(temp %in% NA)]))
		else return(NULL)
	} else {
		options(warn=0)
		stop('Only Entrez Gene or DOLite IDs are supported!')
	}
}

.treeMerge <- function(treeList) {
	if (!is.list(treeList)) stop('Input should be a list! Aborting ...')
	lengthOfCurrentList <- length(treeList) 
	if (lengthOfCurrentList < 2) return(treeList)
	else {
		currentList <- treeList
		while (TRUE) {
			startElmnt <- currentList[[1]]
			remainingElmnt <- lapply(currentList[-1], intersect, startElmnt)
			tempIndex <- which(sapply(remainingElmnt, length) > 0)
			if (length(tempIndex) > 0) {
				tempList <- list(unique(unlist(currentList[c(1,1+tempIndex)])))
				names(tempList) <- paste(names(currentList[c(1,1+tempIndex)]), collapse='***')
				noMerging <- which(sapply(remainingElmnt, length) == 0)
				if (length(noMerging) == 0) {
					currentList <- tempList
					break
				} else {
					currentList <- c(currentList[1+noMerging], tempList)
					lengthOfCurrentList <- length(currentList)
				}
			} else {
				currentList <- c(currentList[-1],currentList[1])
				lengthOfCurrentList <- lengthOfCurrentList - 1
				if (lengthOfCurrentList == 0) break 
			}
		}
		return(currentList)
	}
}

.list2matrix <- function(inputList, verbose=TRUE, removeNA=TRUE) {
	if (!(all(!(is.null(names(inputList)))))) stop('Names of inputlist shoule not contain NULL!')
	if (verbose) print('Removing element containing NULL')
	inputList <- inputList[which(sapply(inputList, length) != 0)]
	column1 <- c()
	if (verbose & removeNA) {
		print('Removing NA')
		removeNA <- function(x) { x <- x[!(x %in% NA)]; return(as.character(unique(x)))}
		tempList <- lapply(inputList, removeNA) 
		tempList <- tempList[which(sapply(tempList, length) != 0)]
	}
	if (length(tempList) > 0) {
		for (i in 1:length(tempList)) {
			column1 <- c(column1, rep(names(tempList)[i], length(tempList[[i]])))
		}
		edgeMatrix <- matrix(column1)
		edgeMatrix <- cbind(edgeMatrix, unlist(tempList))
		edgeMatrix <- edgeMatrix[!(duplicated(edgeMatrix)),]
		if (!is.matrix(edgeMatrix)) edgeMatrix <- matrix(edgeMatrix, ncol=2)
		rownames(edgeMatrix) <- NULL
		return(edgeMatrix)
	} else {
		print('list is empty!')
		return(NULL)
	}
}


.matrix2list <- function(inputMatrix, verbose=TRUE, removeNA=TRUE) {
	if (is.null(inputMatrix)) stop('matrix is NULL!')
	if (verbose & removeNA) {
		print('Removing NA')
		inputMatrix <- inputMatrix[!(is.na(inputMatrix[,1]) | is.na(inputMatrix[,2])) ,]
	}
	if (dim(inputMatrix)[1] > 0) {
		if (!is.matrix(inputMatrix)) {
			inputMatrix <- matrix(inputMatrix, ncol=2)
		}
		newTempList <- as.list(inputMatrix[,2])
		names(newTempList) <- inputMatrix[,1]
		tempNames <- unique(inputMatrix[,1])
		names(tempNames) <- tempNames
		outputList <- lapply(tempNames, function(x, y) {return(unique(c(unlist(y[names(y) %in% x]))))}, newTempList)
		return(outputList)
	} else {
		print('matrix is empty!')
		return(NULL) 
	}

}

