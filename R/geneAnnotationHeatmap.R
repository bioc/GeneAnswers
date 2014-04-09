`geneAnnotationHeatmap` <-
function(annotationList, dataMatrix=NULL, addGeneLabel=TRUE, colorMap=c('#000000', '#FFFFFF'), sortBy='both', standardize.data=TRUE, colorMap.data='GBR', showGeneMax=200, 
			sortBy.data='row', mar=c(1,1,8,6), cex.axis=c(0.8, 0.8), mapType=c('table', 'heatmap'), displayAll=FALSE, symmetry=FALSE, colorBar=FALSE, colorBarLabel=NULL) {
	mapType <- match.arg(mapType) 
	if (colorMap.data[1] == 'GBR') {
		colorMap.data <- rev(RGBColVec(256))
	}	
	allProbe <- unique(unlist(annotationList))
	if (!is.null(dataMatrix)) {
		allProbe <- intersect(allProbe, rownames(dataMatrix))
		# only keep the probes included in the dataMatrix
		annotationList <- lapply(annotationList, function(x) x[x %in% allProbe])
		if (displayAll) dataMatrix <- dataMatrix[rownames(dataMatrix) %in% allProbe,]
		else dataMatrix <- dataMatrix[allProbe,]
	}
	
	if (length(allProbe) > showGeneMax) addGeneLabel <- FALSE
	
	if (!is.null(dataMatrix)) {
		layout(matrix(c(1,2), nrow=1),  widths=c(min(5, max(0.2, 2 * ncol(dataMatrix)/length(annotationList))),1)) 
		#layout(matrix(c(1,2), nrow=1))
		temp <- unique(intersect(unique(unlist(annotationList)), rownames(dataMatrix)))
		newMar <- mar
		newMar[3] <- 2 + ceiling((max(nchar(names(annotationList)))) / 3) 
		newMar[4] <- 1 + ceiling(log2(max(nchar(temp))))
		ord <- .heatmap.mds(dataMatrix, rotate=F, sortBy=sortBy.data, standardize=standardize.data, colorMap=colorMap.data, mar=newMar, cex.axis=cex.axis, 
							addRowLabel=addGeneLabel, mapType='heatmap', symmetry=symmetry, colorBar=colorBar, colorBarLabel=colorBarLabel)
		ord <- ord$row
		if (sortBy != 'none') sortBy <- 'column'
	} else {
		ord <- seq(allProbe)
	}
    
	genes <- unique(unlist(annotationList))
	gene2Annotation <- matrix(0, nrow=length(genes), ncol=length(annotationList))
	rownames(gene2Annotation) <- genes
	colnames(gene2Annotation) <- names(annotationList)
	temp <- sapply(seq(annotationList), function(i) {
		gene.i <- annotationList[[i]]
		gene2Annotation[rownames(gene2Annotation) %in% gene.i,i] <<- 1
	})

	if(!is.null(dataMatrix)) {
		newMar[4] <- 1
		.heatmap.mds(gene2Annotation[ord,], rotate=F, sortBy=sortBy, standardize=F, colorMap=colorMap, mar=newMar, cex.axis=cex.axis, addRowLabel=FALSE, 
				mapType = mapType)
	} else {
		newMar <- c(2,6,(2 + ceiling((max(nchar(names(annotationList)))) / 3)),2)
		if (!addGeneLabel) newMar[2] <- 2
		.heatmap.mds(gene2Annotation[ord,], rotate=F, sortBy=sortBy, standardize=F, colorMap=colorMap, cex.axis=cex.axis, mar=newMar, addRowLabel=addGeneLabel, 
				labelRight=FALSE, mapType = mapType)
	}
   
	if (!is.null(dataMatrix)) layout(matrix(1))
}
