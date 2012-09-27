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
				if ((GOCategory == 'ALL') | (GOCategory == 'BP')) {
					filterGOIDs <- unique(c(filterGOIDs, unlist(lookUp(filterGOIDs, "GO", "BPCHILDREN"))))
					filterGOIDs <- unique(filterGOIDs[!(filterGOIDs %in% NA)])
				}
				if ((GOCategory == 'ALL') | (GOCategory == 'CC')) filterGOIDs <- {
					unique(c(filterGOIDs, unlist(lookUp(filterGOIDs, "GO", "CCCHILDREN"))))
					filterGOIDs <- unique(filterGOIDs[!(filterGOIDs %in% NA)])
				}
				if ((GOCategory == 'ALL') | (GOCategory == 'MF')) filterGOIDs <- {
					unique(c(filterGOIDs, unlist(lookUp(filterGOIDs, "GO", "MFCHILDREN"))))
					filterGOIDs <- unique(filterGOIDs[!(filterGOIDs %in% NA)])
				}
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
	require(igraph0)
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

.drawTable <- function(dataMatrix, mar=c(1,1,5,8), addRowLabel=TRUE, cex.axis=c(1.1, 0.9)) {
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
	cex.axis=c(0.9, 0.9), maxQuantile=0.99, maxVal=NULL, symmetry=FALSE, addRowLabel=TRUE, labelRight=TRUE, rm.unannotate=TRUE, reverseSort=c('none', 'row', 'column', 'both'), 
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
	
	
	if (mapType == 'table') .drawTable(dataMatrix, addRowLabel=addRowLabel, mar=mar, cex.axis=cex.axis)
	
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
			if (labelRight) axis(4, at=1:nrow(dataMatrix), labels=rowLabels, tick=F, las=2, cex.axis=cex.axis[length(cex.axis)])
			else axis(2, at=1:nrow(dataMatrix), labels=rowLabels, tick=F, las=2, cex.axis=cex.axis[length(cex.axis)])	
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
		hyperGTest <- as.data.frame(t(matrix(unlist(testResult), ncol=length(testResult), nrow=length(testResult[[1]]))), stringsAsFactors=FALSE)
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
	require('GO.db')
	temp <- getGOTerm(GOIDs)
	names(temp) <- NULL
	return(unlist(temp)[GOIDs])
}
 
.colorMatch <- function(values, colorMap=NULL, matchMode=c('absolute', 'relative'), zeroColorIndex=NULL, zeroPoint=0) {
	matchMode <- match.arg(matchMode)
	tempColor <- c()

	colorMapping <- function(inputValues, inputColors) {
		tempRange <- max(inputValues)-min(inputValues)
		if (tempRange > 0) scaledValues <- scale(inputValues, scale=(tempRange))
		else scaledValues <- rep(0, length(inputValues))
		outputColors <- inputColors[as.integer((scaledValues - min(scaledValues)) * (length(inputColors) - 1)) + 1]
		names(outputColors) <- as.character(inputValues)
		return(outputColors)
	}

	if (matchMode == 'absolute') {
		if (is.numeric(zeroColorIndex) & (length(zeroColorIndex) == 1)) {
			if (length(unique(values) > 1) | ((zeroColorIndex < length(colorMap)) & (zeroColorIndex > 1))) {
				tempColor <- c(tempColor, colorMapping(c(zeroPoint, values[values >= zeroPoint]), colorMap[zeroColorIndex:length(colorMap)])[-1])
				tempColor <- c(tempColor, colorMapping(c(zeroPoint, values[values < zeroPoint]), colorMap[1:(zeroColorIndex-1)])[-1])
				outColor <- tempColor[as.character(values)]
			} else {
				if ((zeroColorIndex == length(colorMap)) | (zeroColorIndex == 1)) {
					if (unique(values) == zeroPoint) outColor <- colorMap[zeroColorIndex]
					else {
						if (unique(values) > zeroPoint) outColor <- colorMap[length(colorMap)]
						else outColor <- colorMap[1]
					}
				names(outColor) <- as.character(values) 
				} else stop('zeroColorIndex is out of colorMap index! Aborting ...')
			}
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
		stop('Only Entrez Gene or DOLITE IDs are supported!')
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

#.list2matrix <- function(inputList, verbose=TRUE, removeNA=TRUE) {
#	if (!(all(!(is.null(names(inputList)))))) stop('Names of inputlist shoule not contain NULL!')
#	if (verbose) print('Removing element containing NULL')
#	inputList <- inputList[which(sapply(inputList, length) != 0)]
#	column1 <- c()
#	tempList <- inputList
#	if (verbose & removeNA) {
#		print('Removing NA')
#		removeNA <- function(x) { x <- x[!(x %in% NA)]; return(as.character(unique(x)))}
#		tempList <- lapply(inputList, removeNA) 
#		tempList <- tempList[which(sapply(tempList, length) != 0)]
#	}
#	if (length(tempList) > 0) {
#		for (i in 1:length(tempList)) {
#			column1 <- c(column1, rep(names(tempList)[i], length(tempList[[i]])))
#		}
#		edgeMatrix <- matrix(column1)
#		edgeMatrix <- cbind(edgeMatrix, unlist(tempList))
#		edgeMatrix <- edgeMatrix[!(duplicated(edgeMatrix)),]
#		if (!is.matrix(edgeMatrix)) edgeMatrix <- matrix(edgeMatrix, ncol=2)
#		rownames(edgeMatrix) <- NULL
#		return(edgeMatrix)
#	} else {
#		print('list is empty!')
#		return(NULL)
#	}
#}
#
#
#.matrix2list <- function(inputMatrix, verbose=TRUE, removeNA=TRUE) {
#	if (is.null(inputMatrix)) stop('matrix is NULL!')
#	if (verbose & removeNA) {
#		print('Removing NA')
#		inputMatrix <- inputMatrix[!(is.na(inputMatrix[,1]) | is.na(inputMatrix[,2])) ,]
#		if (!is.matrix(inputMatrix)) inputMatrix <- matrix(inputMatrix, ncol=2)
#	}
#	if (dim(inputMatrix)[1] > 0) {
#		if (!is.matrix(inputMatrix)) {
#			inputMatrix <- matrix(inputMatrix, ncol=2)
#		}
#		newTempList <- as.list(inputMatrix[,2])
#		names(newTempList) <- inputMatrix[,1]
#		tempNames <- unique(inputMatrix[,1])
#		names(tempNames) <- tempNames
#		outputList <- lapply(tempNames, function(x, y) {return(unique(c(unlist(y[names(y) %in% x]))))}, newTempList)
#		return(outputList)
#	} else {
#		print('matrix is empty!')
#		return(NULL) 
#	}
#}

.list2matrix <- function(inputList, verbose=TRUE, removeNA=TRUE) {
	if (!(all(!(is.null(names(inputList)))))) stop('Names of the inputlist shoule not contain NULL!')
	if (verbose) print('Removing NULL element ...')
	inputList <- inputList[which(sapply(inputList, length) != 0)]
	
	if (removeNA) {
		if (verbose) print('Removing NA ..., including NA in names of the input list')
		inputList <- inputList[!is.na(names(inputList))]
		inputList <- lapply(inputList, function(x) return(x[!is.na(x)]))
		inputList <- inputList[which(sapply(inputList, length) != 0)]
	}

	if (length(inputList) > 0) {
		if (verbose) print('Converting to a matrix ...')
		inputListLength <- sapply(inputList, length)
		outputMatrix <- matrix(unlist(lapply(names(inputList), function(x, y) return(rep(x,y[x])), inputListLength)))
		names(inputList) <- NULL
		temp <- unlist(inputList)
		outputMatrix <- cbind(outputMatrix, temp)
		if (is.null(names(temp))) {
			rownames(outputMatrix) <- NULL
			outputMatrix <- outputMatrix[!(duplicated(outputMatrix)),]
		} else {
			rownames(outputMatrix) <- names(temp)
			outputMatrix <- outputMatrix[!(duplicated(cbind(outputMatrix, rownames(outputMatrix)))),]
		}
		if (!is.matrix(outputMatrix)) return (t(matrix(outputMatrix)))
		else return(outputMatrix)
	} else {
		print('The input list, probably after removing NA(s), is empty!')
		return(NULL)
	}
}

.matrix2list <- function(inputMatrix, verbose=TRUE, removeNA=TRUE) {
	if (!is.matrix(inputMatrix)) stop('The input is not a matrix. Aborting ...')
	naList <- NULL
	if (removeNA) {
		if (verbose) print('Removing NA ...')
		nonNAIndex <- which(!(is.na(inputMatrix[,1]) | is.na(inputMatrix[,2])))
		if (length(nonNAIndex) == 0) return(NULL)
		if (length(nonNAIndex) == 1) {
			if (verbose) print('Converting to a list ...')
			outputList <- list(inputMatrix[nonNAIndex,2])
			names(outputList) <- inputMatrix[nonNAIndex,1]
			names(outputList[[1]]) <- rownames(inputMatrix)[nonNAIndex]
			return(outputList)
		} else {
			inputMatrix <- inputMatrix[nonNAIndex,]
		}
	} else {
		naIndex <- which(is.na(inputMatrix[,1]))
		if (length(naIndex) > 0) {
			naList <- list(inputMatrix[naIndex, 2])
			names(naList) <- NA
			if (length(naIndex) == 1) names(naList[[1]]) <- rownames(inputMatrix)[naIndex]
		}
	}
	
	if (dim(inputMatrix)[1] > 0) {
		if (verbose) print('Converting to a list ...')
		outputList <- tapply(inputMatrix[,2], inputMatrix[,1], function(i)i, simplify=FALSE)
		listNames <- names(outputList)
		dim(outputList) <- NULL
		names(outputList) <- listNames
		if (!is.null(naList)) outputList <- c(outputList, naList)
		return(outputList) 
	} else {
		if (verbose) print('The input matrix, probably after removing NA(s), is empty!')                     
		return(NULL) 
	}
}

.drawHTMLtable <- function(dataMatrix, outFile, tableName='Table', tableLink= NULL, externalLinkOfTable = NULL, tableCenter=TRUE, catType=c('GO', 'KEGG', 'Entrez', 'DOLITE', 'REACTOME.PATH', 'CABIO.PATH', 'Unknown'), 
						species=c('org.Hs.eg.db', 'org.Rn.eg.db', 'org.Mm.eg.db', 'org.Dm.eg.db'), 
						lastRowLink=FALSE, highlightLastRow=TRUE, matrixOfHeatmap=NULL, topCat=10,  IDCols=1, heatMap=FALSE, reverseOfCluster=FALSE, displayText=TRUE) {
								
	catType <-match.arg(catType)
	species <- match.arg(species)
	
	
	if (lastRowLink & highlightLastRow) stop('Both lastRowLink and highlightLastRow can not be TRUE at the same time!')
	
	if (!(IDCols == 0 | IDCols == 1 | (IDCols == 2 & catType == 'Entrez') | (IDCols == 0 & catType == 'Unknown'))) 
		stop('IDCols must be 0, 1 and 2 (only for Entrez). When catType is Unknown, IDCols must be 0!')
		
	if (is.matrix(dataMatrix) | is.data.frame(dataMatrix) | is.vector(dataMatrix)) {
		dataMatrix <- as.matrix(dataMatrix)
		if ((NA %in% rownames(dataMatrix)) | is.null(rownames(dataMatrix))) print('Warning: NA or NULL might be in rownames of dataMatrix!')
		if ((NA %in% colnames(dataMatrix)) | is.null(colnames(dataMatrix))) print('Warning: NA or NULL might be in colnames of dataMatrix!')
	} else stop('Input is not a valid matrix!')
	if (!(is.null(matrixOfHeatmap))) {
		if (is.matrix(matrixOfHeatmap) | is.data.frame(matrixOfHeatmap)) {
			originalIndexMatrix <- as.matrix(matrixOfHeatmap)
			matrixOfHeatmap <- sqrt(-log10(originalIndexMatrix))
			if ( topCat > 0) matrixOfHeatmap <- apply(matrixOfHeatmap, 2, function(x, top) {if (length(x) > top) x[x < sort(x, decreasing=T)[top]] <- 0; return(x)}, topCat)
		} else stop('matrixOfHeatmap is not a valid matrix!')
	}
	
	if (heatMap & !is.null(matrixOfHeatmap)) {
		conceptCol <- colorRampPalette(c('#77ff00','#ffffff'))
		colorLevel <- 64
		maxColor <- max(sqrt(-log10(originalIndexMatrix)))
		minColor <- min(sqrt(-log10(originalIndexMatrix)))
		contentsOfColorlevel <- conceptCol(colorLevel)
		if (!reverseOfCluster) contentsOfColorlevel <- rev(contentsOfColorlevel)
	}
		 
	HTwrap <- function(x, tag = "TD", scripts=NULL) {
        if (is.null(scripts)) return(paste("<", tag, ">", x, "</", tag, ">", sep = ""))
		else return(paste("<", tag, ' ', scripts, ">", x, "</", tag, ">", sep = ""))
    }
    
    cat(paste('<H2 align=center><font face="courier" size="2">', sep=''), file = outFile, sep = "")
    if (is.null(tableLink))	{
		if (is.null(externalLinkOfTable)) cat(paste(tableName, "</font></H2>", sep=''), file = outFile, sep = "\n")
		else cat(paste('<a href="', externalLinkOfTable, '">',tableName, "</a></font></H2>", sep=''), file = outFile, sep = "\n")
	} else {
		if (is.null(externalLinkOfTable)) cat(paste('<a name="', tableLink, '">', tableName, "</a></font></H2>", sep=''), file = outFile, sep = "\n")
		else cat(paste('<a name="', tableLink, '" href="', externalLinkOfTable, '">', tableName, "</a></font></H2>", sep=''), file = outFile, sep = "\n")
	}
	if (tableCenter) cat("<center> \n", file = outFile)
    cat('<TABLE BORDER=1, style="font-family:courier;text-align:center;font-size:12px">', file = outFile, sep = "\n")
	
	if (IDCols == 0) {
		if (catType == 'Unknown') firstLine <- HTwrap('Customized IDs')
		else firstLine <- HTwrap(paste(catType, ' Terms/Names', sep=''))
	} else {
		if (IDCols == 1) {
			if (catType == 'Entrez')  firstLine <- paste(HTwrap(paste('Gene Symbols', sep='')), HTwrap(paste(catType, ' IDs', sep='')), sep='')
			else {
				if (catType == 'Unknown') firstLine <- HTwrap(paste(catType, ' Terms', sep=''))
				else firstLine <- paste(HTwrap(HTwrap(paste(catType, ' Terms', sep=''), tag='A', scripts='HREF="#h-2"')), HTwrap(paste(catType, ' IDs', sep='')),sep='')
			}
		}
		else firstLine <- paste(HTwrap(paste('Gene Symbols', sep='')), HTwrap(paste(catType, ' IDs', sep='')), HTwrap(paste('Gene Names', sep='')),sep='')
	}
	
	
	if (catType == 'Entrez') {
		firstLine <- paste(firstLine, paste(HTwrap(colnames(dataMatrix)), collapse=''), sep='')
		indexM <- NULL
	} else {
		firstLine <- paste(firstLine, paste(HTwrap(HTwrap(colnames(dataMatrix), tag='A', scripts=paste('HREF="#h-', 2+c(1:dim(dataMatrix)[2]), '"', sep=''))), sep='', collapse=''), sep='')
		indexM <- which((is.numeric(dataMatrix) & (dataMatrix > 0)) | (is.character(dataMatrix) & (dataMatrix != "0 (1)")), arr.ind = TRUE)
	}
	
	cat(firstLine, '\n', file = outFile, sep="")
	
	if (highlightLastRow) {
		rownames(dataMatrix)[dim(dataMatrix)[1]] <- paste('<b><i>', rownames(dataMatrix)[dim(dataMatrix)[1]], '</i></b>', sep='')
	}
	
	rowIDs <- rownames(dataMatrix)
	if (lastRowLink) {
		linkRowIDs <- length(rowIDs)
	} else {
		linkRowIDs <- length(rowIDs) - 1
	}
	
	geneMapping <- function(geneIDs, data, info=c('SYMBOL', 'GENENAME')) {
		info <- match.arg(info)
        temp <- unlist(lookUp(geneIDs, data, info))
		temp[temp %in% NA] <- 'Discontinued'
		return(temp)
	}
	if (IDCols == 0) {
		hyperLinkPrefix <- switch(catType,
			'GO'=paste('<a href=http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=', rowIDs[1:linkRowIDs], '>', getCategoryTerms(rowIDs[1:linkRowIDs], catType=catType, missing='keep'), '</a>',sep=''),
			'KEGG'=paste('<a href=http://www.genome.jp/dbget-bin/www_bget?ko', rowIDs[1:linkRowIDs], '>', getCategoryTerms(rowIDs[1:linkRowIDs], catType=catType, missing='keep'), '</a>',sep=''),
			'DOLITE'=paste('<a>', getCategoryTerms(rowIDs[1:linkRowIDs], catType=catType, missing='keep'), '</a>',sep=''),
			'REACTOME.PATH'=paste('<a href=http://www.reactome.org/cgi-bin/eventbrowser?DB=gk_current&ID=', rowIDs[1:linkRowIDs], '&>', getCategoryTerms(rowIDs[1:linkRowIDs], catType=catType, missing='keep'), '</a>',sep=''),
			'CABIO.PATH'=paste('<a href=', unlist(lapply(rowIDs[1:linkRowIDs], .caBIOLink)), '>', getCategoryTerms(rowIDs[1:linkRowIDs], catType=catType, missing='keep'), '</a>', sep=''),
			'Entrez'=paste('<a href=http://www.ncbi.nlm.nih.gov/sites/entrez?Db=gene&Cmd=ShowDetailView&TermToSearch=', rowIDs[1:linkRowIDs], '>', geneMapping(rowIDs[1:linkRowIDs], species, info='SYMBOL'), '</a>',sep=''),
			'Unknown'=rowIDs[1:linkRowIDs])
	} else {
		hyperLinkPrefix <- switch(catType,
			'GO'=paste(getCategoryTerms(rowIDs[1:linkRowIDs], catType=catType, missing='keep'), '</TD><TD><a href=http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=', rowIDs[1:linkRowIDs], '>',  rowIDs[1:linkRowIDs], '</a>',sep=''),
			'KEGG'=paste(getCategoryTerms(rowIDs[1:linkRowIDs], catType=catType, missing='keep'), '</TD><TD><a href=http://www.genome.jp/dbget-bin/www_bget?ko', rowIDs[1:linkRowIDs], '>', rowIDs[1:linkRowIDs], '</a>',sep=''),
			'DOLITE'=paste(getCategoryTerms(rowIDs[1:linkRowIDs], catType=catType, missing='keep'), '</TD><TD><a>', rowIDs[1:linkRowIDs], '</a>',sep=''),
			'REACTOME.PATH'=paste(getCategoryTerms(rowIDs[1:linkRowIDs], catType=catType, missing='keep'), '</TD><TD><a href=http://www.reactome.org/cgi-bin/eventbrowser?DB=gk_current&ID=', rowIDs[1:linkRowIDs], '&>', rowIDs[1:linkRowIDs], '</a>',sep=''),
			'CABIO.PATH'=paste(getCategoryTerms(rowIDs[1:linkRowIDs], catType=catType, missing='keep'), '</TD><TD><a href=', unlist(lapply(rowIDs[1:linkRowIDs], .caBIOLink)), '>', rowIDs[1:linkRowIDs], '</a>',sep=''), 
			'Entrez'=paste(geneMapping(rowIDs[1:linkRowIDs], species, info='SYMBOL'), '</TD><TD><a href=http://www.ncbi.nlm.nih.gov/sites/entrez?Db=gene&Cmd=ShowDetailView&TermToSearch=', rowIDs[1:linkRowIDs], '>', rowIDs[1:linkRowIDs], '</a>',sep=''),
			'Unknown'=rowIDs[1:linkRowIDs])
		if (IDCols > 1) hyperLinkPrefix <- paste(hyperLinkPrefix, '</TD><TD>', geneMapping(rowIDs[1:linkRowIDs], species, info='GENENAME'), sep='')
	}
	
   	
    if (!lastRowLink) {
		if (IDCols == 0) hyperLinkPrefix <- c(hyperLinkPrefix, rowIDs[length(rowIDs)])
		else {
			if (IDCols == 1) { 
				if (catType == 'Unknown') hyperLinkPrefix <- c(hyperLinkPrefix, paste(rowIDs[length(rowIDs)], sep=''))
				else hyperLinkPrefix <- c(hyperLinkPrefix, paste(rowIDs[length(rowIDs)], '</TD><TD>',sep=''))
			}
			else hyperLinkPrefix <- c(hyperLinkPrefix, paste(rowIDs[length(rowIDs)], '</TD><TD></TD><TD>',sep=''))
		}
	}
    
	

	if (topCat > 0 & !is.null(matrixOfHeatmap)) {
		matrixOfHeatmapIndex <- which(matrixOfHeatmap > 0, arr.ind=TRUE)
	} else matrixOfHeatmapIndex <- NULL
	
    if (is.null(matrixOfHeatmap)) sigIndex <- NULL
	else  sigIndex <- which(originalIndexMatrix < 1, arr.ind=TRUE)
	
    if (!displayText) dataMatrix[1:(dim(dataMatrix)[1]-1),] <- '  '

	for (i in 1:length(rowIDs)) {
		textLine <- HTwrap(hyperLinkPrefix[i])
		for (j in 1:length(colnames(dataMatrix))) {
			if (dataMatrix[i,j] == '') dataMatrix[i,j] <- 'NA'
			if (i == dim(dataMatrix)[1] & highlightLastRow)  dataMatrix[i,j] <- HTwrap(HTwrap(dataMatrix[i,j], tag='b'), tag='i') 
			if (is.null(indexM)) textLine <- paste(textLine, HTwrap(dataMatrix[i,j]), sep='')
			else {
				if ((i %in% indexM[indexM[,2] == j,1])){
					k <- dim(dataMatrix)[2] + match(i, indexM[indexM[,2] == j,1]) + match(j, indexM[,2])
					bgColor <- NULL
					# if heatmap is on, bgcolor <- contentsOfColorlevel[ceiling(colorlevel * sqrt(-log10(originalIndexMatrix)[i,j])/(sqrt(-log10(max)) - sqrt(-log10(min)) + 1))]
					if (heatMap) {
						if (!is.null(sigIndex) & (i %in% sigIndex[sigIndex[,2] == j,1])) {
							dataMatrix[i,j] <- HTwrap(dataMatrix[i,j], tag='b')
							bgColor <- paste('bgcolor="', contentsOfColorlevel[ceiling((colorLevel-1) * sqrt(-log10(originalIndexMatrix)[i,j])/(maxColor - minColor)) + 1], '"', sep='')
						}
					} else {
						if (!is.null(sigIndex) & (i %in% sigIndex[sigIndex[,2] == j,1])) dataMatrix[i,j] <- HTwrap(dataMatrix[i,j], tag='b')
						if (!is.null(matrixOfHeatmapIndex) & (i %in% matrixOfHeatmapIndex[matrixOfHeatmapIndex[,2] == j,1])) bgColor <- 'bgcolor="#ccff00"'
					}
					
					textLine <- paste(textLine, HTwrap(HTwrap(dataMatrix[i,j], tag='A', scripts=paste('href="#h-', (k+1), '"', sep='')), scripts=bgColor), sep='')
				}else textLine <- paste(textLine, HTwrap(dataMatrix[i,j]), sep='')
			}			
		}
		cat(HTwrap(textLine, tag = 'TR'), '\n', file = outFile, sep="")
	}
	
	cat("</TABLE>", sep = "\n", file = outFile)
    if (tableCenter) cat("</center> \n", file = outFile)
}

.name2lib <- function(species=c('anopheles', 'arabidopsis', 'bovine', 'worm', 'canine', 'fly', 'zebrafish', 'ecolistraink12', 'ecolistrainsakai', 'chicken', 'human', 'mouse', 'rhesus', 'malaria', 'chimp', 
								'rat', 'yeast', 'pig', 'xenopus')) {
	return(switch(tolower(species), 
		'anopheles'='org.Ag.eg.db', 'arabidopsis'='org.At.tair.db', 'bovine'='org.Bt.eg.db', 'worm'='org.Ce.eg.db', 'canine'='org.Cf.eg.db', 'fly'='org.Dm.eg.db', 'zebrafish'='org.Dr.eg.db',
		'ecolistraink12'='org.EcK12.eg.db', 'ecolistrainsakai'='org.EcSakai.eg.db', 'chicken'='org.Gg.eg.db', 'human'='org.Hs.eg.db', 'mouse'='org.Mm.eg.db', 'rhesus'='org.Mmu.eg.db', 
		'malaria'='org.Pf.plasmo.db', 'chimp'='org.Pt.eg.db', 'rat'='org.Rn.eg.db', 'yeast'='org.Sc.sgd.db', 'pig'='org.Ss.eg.db', 'xenopus'='org.Xl.eg.db', NULL))
}

#kernal xml query function
.getcaBIOIDInfo <- function(xmlLink, IDinfo=c('indirect', 'direct', 'path', 'geneNumber')) {
	IDinfo <- match.arg(IDinfo)
	require(XML)
	root <- xmlChildren(xmlRoot(try(xmlTreeParse(file=URLencode(xmlLink), isURL=TRUE))))
	if ('queryResponse' %in% names(root)) {
		level1 <- xmlChildren(root[['queryResponse']])
		counts <- as.numeric(xmlValue(level1[['recordCounter']]))
		if (IDinfo == 'geneNumber') return(counts)
		if ( counts > 0) {
			getID <- function(x, IDinfo) {
				temp <- xmlChildren(x)
				content <- lapply(temp, xmlAttrs)
				names(content) <- lapply(temp, xmlValue)
				if (IDinfo == 'indirect') {
					catch <- function(y) { if (y['name'] == 'id') {names(y) <- NULL; return(y)}}
					return(names(unlist(sapply(content, catch))))
				}
				if (IDinfo == 'direct') {
					tempVector <- unlist(content)
					if (names(tempVector[tempVector == 'sourceType']) == 'Entrez gene.name') {
						return(sub(".name", "", names(tempVector[tempVector == "crossReferenceId"])))
					}
				}
				if (IDinfo == 'path') {
					tempVector <- unlist(content)
					tempList <- lapply(c("displayValue", "description", "name", "source"), function(x, y) {return(sub(".name", "", names(y[y==x])))}, tempVector)
					names(tempList) <- c("displayValue", "description", "name", "source")
	   				return(unlist(tempList))
				}
			}
			result <- unlist(sapply(level1[which(names(level1) == "class")], getID, IDinfo))
			if (is.matrix(result)) colnames(result) <- NULL
			else names(result) <- NULL
			return(result)
		}
	}
}

# retrieve caBIO pathway info
.getcaBIOPATHInfo <- function(caBIOPATHIDs) {
	.getPathInfo <- function (x, ...) {
		if (is.na(as.numeric(x))) { 
			return(NA)
		} else {
			pathInfo <- unlist(as.list(.getcaBIOIDInfo(paste('http://cabioapi.nci.nih.gov/cabio43/GetXML?query=Pathway[@id=',x, ']', sep=''), ...)[,1]))
			if (is.null(pathInfo)) return(NA)
			else return(pathInfo)
		}
	}
	temp <- lapply(caBIOPATHIDs, .getPathInfo, IDinfo='path')
	names(temp) <- caBIOPATHIDs
	return(temp)
}

.getcaBIOPATHs <- function(caBIOIDs) {
	print('Retrieving pathways in caBIO ...')
	temp <- lapply(paste('http://cabioapi.nci.nih.gov/cabio43/GetXML?query=Pathway&Gene[@id=', unlist(caBIOIDs), ']',sep=''), .getcaBIOIDInfo)
	names(temp) <- unlist(caBIOIDs)
	return(temp)
}

#retrieve all caBIO genes based on the given caBIO pathway IDs 
#"http://cabioapi.nci.nih.gov/cabio43/GetXML?query=Gene&Pathway[@id=95]"
.getcaBIOPATHGenes <- function(caBIOPATHIDs) {
	print('Retrieving pathways genes in caBIO ...')
	temp <- lapply(paste('http://cabioapi.nci.nih.gov/cabio43/GetXML?query=Gene&Pathway[@id=', unlist(caBIOPATHIDs), ']',sep=''), .getcaBIOIDInfo)
	names(temp) <- unlist(caBIOPATHIDs)
	return(temp)
}


#map all caBIO Gene IDs in a list to entrez IDs
.mapListcaBIOIDs2entrez <- function(caBIOList, hide=TRUE) {
	allMapping <- caBIO2entrez(unique(unlist(caBIOList)))
	mappingVector <- function(x, y, hide=hide) {
		temp <- unique(unlist(y[names(y) %in% x]))
		if (hide) names(temp) <- NULL
		return(temp)
	}
	return(lapply(caBIOList, mappingVector, allMapping, hide=hide))
}

#generate pathway link based on the given caBIO pathway ID
.caBIOLink <- function(caBIOPATHID){
	temp <- .getcaBIOPATHInfo(caBIOPATHID)[[1]]
	if (all(c("displayValue", "description", "name", "source") %in% names(temp))) {
		return(switch(temp['source'],
				'NCI-Nature Curated'=paste('http://pid.nci.nih.gov/search/pathway_landing.shtml?pathway_id=', substr(temp['name'], 3, nchar(temp['name'])), '&pathway_name=', 
											gsub(' ', '%20', temp['displayValue']), '&source=NCI-Nature%20curated&what=graphic&jpg=on', sep=''),
				'Reactome'=paste('http://pid.nci.nih.gov/search/pathway_landing.shtml?source=Reactome%20Imported&what=graphic&jpg=on&pathway_id=', substr(temp['name'], 3, nchar(temp['name'])), sep=''),
				'BioCarta'=paste('http://cgap.nci.nih.gov/Pathways/BioCarta/', temp['name'], sep='')))
	}
}

#generate concepts relationship for given concepts.
.catsCluster <- function(dataMatrix, gAL, clusterMethod=c("ward", "single", "complete", "average", "mcquitty", "median", "centroid"), catTerm=TRUE, 
						catType=c('GO', 'KEGG', 'DOLITE', 'REACTOME.PATH', 'CABIO.PATH', 'Unknown'), ...) {
	clusterMethod <- match.arg(clusterMethod)
	catType <-match.arg(catType)
	if (!is.matrix(dataMatrix)) dataMatrix <- matrix()
	catsGenes <- lapply(rownames(dataMatrix), function(x, y) {return(unique(unlist(lapply(lapply(y, getGenesInCategory), function(m,n) {return(m[names(m) %in% n])}, x))))}, gAL)
	names(catsGenes) <- rownames(dataMatrix) 
	allGenes <- unique(unlist(catsGenes))
   #indexCats <- as.character(c(1:length(catsGenes)))
	#names(indexCats) <- rownames(dataMatrix)
	if (catTerm & (catType != 'Unknown')) {
		names(catsGenes) <- getCategoryTerms(names(catsGenes), catType, ...)
	}
	clusterM <- matrix(0, nrow=length(allGenes), ncol=length(catsGenes), dimnames=list(allGenes, names(catsGenes))) 
	for (i in 1:dim(clusterM)[2]) {
		clusterM[catsGenes[[i]],i] <- 1
	}
	hc <- hclust(dist(t(clusterM)), method=clusterMethod)
	#print(ceiling(log2(max(nchar(colnames(clusterM))))))
    op <- par(mar=c(1,1,1,40))
	plot(as.dendrogram(hc), horiz=TRUE, xlab='Concepts', ylab=NULL, main=NULL, axes=FALSE, ann=FALSE)
	par(op)
	return(invisible(clusterM))
}

# graphAttr should be a 2 element list, each one is a dataframe, represents vertex or edge attributes, repectively
.convertCytoscapeWeb <- function(graphAttr, htmlName=NULL,  fileSuffix=1, verbose=TRUE, destination=NULL, bgColor='#ffffff') {
	if (is.null(htmlName)) cytoDir <- 'cytoscapeWebFiles'
	else cytoDir <- paste(htmlName, 'cytoscapeWebFiles', sep='.')
	if (dir.create(cytoDir, showWarnings=FALSE) | file.exists(cytoDir)) {
		setwd(cytoDir)
		if (verbose) print(paste('Generating ', paste('edges_', fileSuffix, '.txt', sep=''), sep=''))
		write.table(graphAttr[['edge.attributes']][,c('NODES1', 'NODES2')], file=paste('edges_', fileSuffix, '.txt', sep=''), quote=FALSE, row.names=FALSE, sep='\t')
		if (verbose) print(paste('Generating ', paste('edgesAttr_', fileSuffix, '.txt', sep=''), sep=''))
		edgeAttrM <- cbind(paste(graphAttr[['edge.attributes']][,'NODES1'], '()', graphAttr[['edge.attributes']][,'NODES2']), graphAttr[['edge.attributes']][, c(3:dim(graphAttr[['edge.attributes']])[2])],stringsAsFactors = FALSE)
		colnames(edgeAttrM)[1] <- c('EDGES')
		write.table(edgeAttrM, file=paste('edgesAttr_', fileSuffix, '.txt', sep=''), quote=FALSE, row.names=FALSE, sep='\t')
		if (verbose) print(paste('Generating ', paste('verticesAttr_', fileSuffix, '.txt', sep=''), sep=''))
		write.table(graphAttr[['vertex.attributes']], file=paste('verticesAttr_', fileSuffix, '.txt', sep=''), quote=FALSE, row.names=FALSE, sep='\t')
		setwd('..')
		
		if (file.exists('NetworkTransformer.jar')) {
			originalWD <- getwd()
			setwd(cytoDir)
			system(paste('java -jar ../NetworkTransformer.jar', paste('edges_', fileSuffix, '.txt', sep=''), paste('edgesAttr_', fileSuffix, '.txt', sep=''), paste('verticesAttr_', fileSuffix, '.txt', sep=''), fileSuffix, substr(bgColor, 2, 7)))
			if(!is.null(destination)) {
				if (file.copy(paste('edges_', fileSuffix, '.js', sep=''), destination)) {
					file.remove(paste('edges_', fileSuffix, '.js', sep=''))
				} else {
					setwd(originalWD)
					stop('Javascript file ', paste('edges_', fileSuffix, '.js', sep=''), ' can not be moved! Working directory is reset to ', originalWD, ' , Aborting ...')
				}
			}
			setwd('..')
			return(invisible(paste('edges_', fileSuffix, '.js', sep='')))
		} else {
			stop('Necessary java program is missing, Aborting ...')
		}
	} else stop('Failure to create ', cytoDir, ' subdirectory! Aborting ...')
}
