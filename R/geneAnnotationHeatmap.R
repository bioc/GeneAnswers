`geneAnnotationHeatmap` <-
function(annotationList, dataMatrix=NULL, addGeneLabel=TRUE, colorMap=c('#000000', '#FFFFFF'), sortBy='both', standardize.data=TRUE, colorMap.data='GBR', sortBy.data='row', mar=c(1,1,8,6), cex.axis=c(0.8, 0.8), mapType=c('table', 'heatmap'), ... ) {
	mapType <- match.arg(mapType)
	if (colorMap.data == 'GBR') {
		require(Heatplus)
		colorMap.data <- rev(RGBColVec(256))
	} else {
		require(RColorBrewer)
	}
	require(MASS)
	
	allProbe <- unique(unlist(annotationList))
	if (!is.null(dataMatrix)) {
		allProbe <- intersect(allProbe, rownames(dataMatrix))
		# only keep the probes included in the dataMatrix
		annotationList <- lapply(annotationList, function(x) x[x %in% allProbe])
		dataMatrix <- dataMatrix[allProbe,]
	}
	
	if (length(allProbe) > 300) addGeneLabel <- F
	
	if (!is.null(dataMatrix)) {
		layout(matrix(c(1,2), nrow=1),  widths=c(min(5, max(0.2, 1.2 * ncol(dataMatrix)/length(annotationList))),1)) 
		#layout(matrix(c(1,2), nrow=1))
		temp <- unique(intersect(unique(unlist(annotationList)), rownames(dataMatrix)))
		newMar <- mar
		newMar[3] <- 2 + ceiling((max(nchar(names(annotationList)))) / 3) 
		newMar[4] <- 1 + ceiling(log2(max(nchar(temp))))
		ord <- .heatmap.mds(dataMatrix, rotate=F, sortBy=sortBy.data, standardize=standardize.data, colorMap=colorMap.data, mar=newMar, cex.axis=cex.axis, addRowLabel=addGeneLabel, mapType='heatmap', ...)
		ord <- ord$row
		sortBy <- 'column'
	} else {
		ord <- seq(allProbe)
	}
	gene2Annotation <- matrix(0, nrow=length(allProbe), ncol=length(annotationList))
	rownames(gene2Annotation) <- allProbe
	colnames(gene2Annotation) <- names(annotationList)
	temp <- sapply(seq(annotationList), function(i) {
		gene.i <- annotationList[[i]]
		gene2Annotation[allProbe %in% gene.i,i] <<- 1
	})

	newMar[4] <- 1
	.heatmap.mds(gene2Annotation[ord,], rotate=F, sortBy=sortBy, standardize=F, colorMap=colorMap, mar=newMar, cex.axis=cex.axis, addRowLabel=FALSE, mapType = mapType, ...)
   
	if (!is.null(dataMatrix)) layout(matrix(1))
}

