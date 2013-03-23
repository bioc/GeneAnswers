`drawTable` <-
function(dataMatrix, topCat=10, heatMap=TRUE, matrixOfHeatmap=NULL, clusterTable=c('geneNum', 'pvalue', NULL), methodOfCluster=c('mds', 'sort'), mar=c(1,5,5,8), addRowLabel=TRUE, cex.axis=c(1.1, 0.9), 
		reverseOfCluster=FALSE, xGridLine=FALSE,  colorBar=TRUE, newWindow=TRUE, endOfColBar=c('> 0.01', 'Minimum of p values'), heatMapColor=c('#00ff00','#ffffff'), canvasWidth=NULL, canvasHeight=NULL, ...) {
	if (!is.null(clusterTable))  {
		clusterTable <- match.arg(clusterTable)
		methodOfCluster <- match.arg(methodOfCluster)
	}
	if (is.matrix(dataMatrix) | is.data.frame(dataMatrix)) {
		if (is.data.frame(dataMatrix)) dataMatrix <- as.matrix(dataMatrix)
	} else stop('Input is not a valid matrix!')
	if (!(is.null(matrixOfHeatmap))) {
		if (is.matrix(matrixOfHeatmap) | is.data.frame(matrixOfHeatmap)) {
			originalIndexMatrix <- as.matrix(matrixOfHeatmap)
			matrixOfHeatmap <- sqrt(-log10(originalIndexMatrix))
		} else stop('matrixOfHeatmap is not a valid matrix!')
	}
	
	tableTitle <- c("Top Categories Table")
	
	if (!is.null(clusterTable) & !is.null(matrixOfHeatmap)) {
		if (clusterTable == 'geneNum') {
			if (is.numeric(dataMatrix)) {
				baseM <- dataMatrix[1:(dim(dataMatrix)[1]-1),]
				tableTitle <- paste(tableTitle, "Culstered by Gene Numbers")
			}
			else stop('The input Matrix is not numeric!')
		} else {
			baseM <- matrixOfHeatmap
			tableTitle <- paste(tableTitle, "Culstered by p Values") 
		}
		
		if (methodOfCluster == 'mds') {
			distBaseM <- as.matrix(dist(baseM, method = "euclidean", upper=TRUE))
			distBaseM[distBaseM < 1e-6] <- 1e-6
			diag(distBaseM) <- 0
			clusterResult <- isoMDS(distBaseM, k=1)$points
			clusterOrder <- rownames(clusterResult)[order(clusterResult, decreasing=reverseOfCluster)]
			originalIndexMatrix <- originalIndexMatrix[clusterOrder,]
			matrixOfHeatmap <- matrixOfHeatmap[clusterOrder,]
			tempM <- dataMatrix[c(clusterOrder, rownames(dataMatrix)[dim(dataMatrix)[1]]),]
			dataMatrix <- tempM
		} else {
			clusterOrder <- sort(rownames(baseM), ...)
			tempM <- dataMatrix[c(clusterOrder, rownames(dataMatrix)[dim(dataMatrix)[1]]),]
			dataMatrix <- tempM
		}
			
	}

	if (heatMap) {
		conceptCol <- colorRampPalette(heatMapColor)
		colorLevel <- 32
	}
	
	if (is.null(rownames(dataMatrix)) | is.null(colnames(dataMatrix))) stop('rownames and/or colnames of input matrix are missing!')
	
	if (newWindow) {
		if (is.null(canvasWidth)) canvasWidth <- 2*(dim(dataMatrix)[2])*ceiling(max(nchar(as.character(dataMatrix)))/5)
		if (is.null(canvasHeight)) canvasHeight <- dim(dataMatrix)[1]/3
		x11(width=canvasWidth, height=canvasHeight, title=tableTitle)
	}
	oldsetting <- par('mar'=mar)
	if (is.null(matrixOfHeatmap)) image(x=1:ncol(dataMatrix), y=1:nrow(dataMatrix), -t(matrix(0, nrow=nrow(dataMatrix), ncol=ncol(dataMatrix)) ), col=c('white','white'), axes=F, xlab='', ylab='', )
	else {
		matrixOfHeatmap <- rbind(matrixOfHeatmap, total=rep(0, length=dim(matrixOfHeatmap)[2]))
 		if (heatMap) {
 			image(1:ncol(matrixOfHeatmap), 1:nrow(matrixOfHeatmap), t(apply(matrixOfHeatmap,2,rev)), col=(rev(conceptCol(colorLevel))), axes=F, xlab='', ylab='')
 			if (topCat > 0) matrixOfHeatmap <- apply(matrixOfHeatmap, 2, function(x, top) {x[x < sort(x, decreasing=T)[top]] <- 0; return(x)}, topCat)
 		} 
 		else {
 			if (topCat > 0) {
				matrixOfHeatmap <- apply(matrixOfHeatmap, 2, function(x, top) {if (length(grep(0, sort(x, decreasing=T))) > 0) top <- min(which(sort(x, decreasing=T) == 0)[1] - 1, top); 
					x[x < sort(x, decreasing=T)[top]] <- 0; x[x >= sort(x, decreasing=T)[top]] <- 1;return(x)}, topCat)
 				colorMap = c('#d0ffd0','#ffffff')
 				image(1:ncol(matrixOfHeatmap), 1:nrow(matrixOfHeatmap), t(apply(matrixOfHeatmap,2,rev)), col=rev(colorMap), axes=F, xlab='', ylab='', )
		   	} else {
			    image(x=1:ncol(dataMatrix), y=1:nrow(dataMatrix), -t(matrix(0, nrow=nrow(dataMatrix), ncol=ncol(dataMatrix)) ), col=c('white','white'), axes=F, xlab='', ylab='', )
			}
 		}
 	}
	if(xGridLine) abline(h=c(0:nrow(dataMatrix))+0.5, col='#bbbbbb', lty = "solid", lwd=2)  
	if (addRowLabel) axis(2,at=1:nrow(dataMatrix), labels=rev(rownames(dataMatrix)), tick=F,las=2, cex.axis=cex.axis[length(cex.axis)])
	axis(3,at=1:ncol(dataMatrix), labels=colnames(dataMatrix), tick=F,las=1, cex.axis=cex.axis[1])
	
	if (topCat == 0) matrixOfHeatmap <- NULL
	if (!is.null(matrixOfHeatmap)) {
		matrixOfHeatmapIndex <- which(matrixOfHeatmap > 0, arr.ind=TRUE)
		text(matrixOfHeatmapIndex[,2], (dim(dataMatrix)[1] - matrixOfHeatmapIndex[,1] + 1), label=dataMatrix[matrixOfHeatmapIndex], col='red', cex=3/log(nrow(dataMatrix)), font=4)
		if (topCat <= (nrow(matrixOfHeatmap) - 2)) {
			matrixOfHeatmapIndex <- which(matrixOfHeatmap[c(1:(dim(matrixOfHeatmap)[1]-1)),] == 0, arr.ind=TRUE)
	  		text(matrixOfHeatmapIndex[,2], (dim(dataMatrix)[1] - matrixOfHeatmapIndex[,1] + 1), label=dataMatrix[matrixOfHeatmapIndex], col='black', cex=3/log(nrow(dataMatrix)), font=1)
		}
	} else {
		for (i in 1:dim(dataMatrix)[2]) {
			text(rep(i, length=dim(dataMatrix)[2]), c(2:dim(dataMatrix)[1]), label=rev(dataMatrix[1:(dim(dataMatrix)[1]-1),i]), col='black', cex=3/log(nrow(dataMatrix)))
		}
	}
	text(c(1:dim(dataMatrix)[2]), rep(1, length=dim(dataMatrix)[2]), label= dataMatrix[dim(dataMatrix)[1],], col='black', cex=3/log(nrow(dataMatrix)), font=4)
	par(oldsetting)

	if(colorBar & heatMap & !is.null(matrixOfHeatmap)) {
		x11(width=2.5, title="Color Bar") 
		oldsetting <- par('mar'=c(3,2,3,10))
		image(1, c(1:colorLevel), matrix(c(1:colorLevel), 1, colorLevel), col=(rev(conceptCol(colorLevel))), axes=F, xlab='', ylab='')
		par(las = 1)
		ruler <- rep('', length=colorLevel)
		minPvalue <- min(originalIndexMatrix)
		if (!is.null(endOfColBar)) {
			ruler[1] <- endOfColBar[1]
			if (length(endOfColBar) > 1)  ruler[colorLevel] <- endOfColBar[2]
			else ruler[colorLevel] <- as.character(signif(minPvalue, digits=3))
		} else {
			ruler[1] <- "> 0.01"
			ruler[colorLevel] <- as.character(signif(minPvalue, digits=3))
		}
		ruler[colorLevel/2] <- as.character(signif(10^(-(sqrt(-log10(minPvalue))/2)^2), digits=3))
		ruler[colorLevel/4] <- as.character(signif(10^(-(sqrt(-log10(minPvalue))/4)^2), digits=3))
		ruler[3 * colorLevel/4] <- as.character(signif(10^(-(3 * sqrt(-log10(minPvalue))/4)^2), digits=3))
		axis(4, at = c(1:colorLevel), labels = ruler, tick=FALSE, cex.axis=3/log(nrow(dataMatrix))) 
        box()
        las <- 0
		par(oldsetting)
		dev.set(dev.prev())
		#x11.options(reset = TRUE)
	}
	#x11.options(reset = TRUE)
}