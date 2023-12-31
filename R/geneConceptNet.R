`geneConceptNet` <-
function(inputList, lengthOfRoots=NULL, inputValue=NULL, centroidSize='geneNum', output=c('fixed','interactive', 'none'), colorMap=NULL, bgColor='#ffffff', matchMode=c('absolute', 'relative'), zeroColorIndex=NULL, verbose=FALSE, symmetry=TRUE) {
	output <- match.arg(output)
	matchMode <- match.arg(matchMode)
	
	if (is.null(lengthOfRoots)) lengthOfRoots <- length(inputList)
	g <- .list2graph(inputList, verbose=verbose)
#	g <- temp[[1]]
#	nodes <- temp[[2]]
#	vertexLabels <- names(nodes)
#	vertexLabels[1:length(inputList)] <- names(inputList) 
	E(g)$color <- "#9999ff"
	E(g)$width <- 3
	V(g)$size <- 8
	V(g)$color <- "#dddddd"
	V(g)$label <- V(g)$name
	if (!is.null(inputValue)) {
		if(all(!is.na(as.numeric(inputValue))) & (all(unique(unlist(inputList)) %in% names(inputValue)))) {
			if (is.null(colorMap)) {
				colorLevel <- 256
#				if ((min(as.numeric(inputValue)) == 0) & (max(as.numeric(inputValue)) == 0))
				if (all(as.numeric(inputValue) == 0)) {
					matchMode='relative'
					colorMap=rep(V(g)$color, length(inputValue))
					zeroColorIndex <- NULL
				} else {
					if (all(as.numeric(inputValue) > 0) | all(as.numeric(inputValue) < 0)) {
						if (all(as.numeric(inputValue) > 0)) {
							zeroColorIndex <- 1
							conceptCol <- colorRampPalette(c(bgColor,'#ff3333'))
						} else {
							zeroColorIndex <- colorLevel
							conceptCol <- colorRampPalette(c('#3333ff',bgColor))
						}
						colorMap <- conceptCol(colorLevel)
						symmetry <- FALSE
					} else {
						#zeroColorIndex <- 1 + ceiling(abs(min(as.numeric(inputValue))) * (colorLevel-1) / (max(as.numeric(inputValue))-min(as.numeric(inputValue))))
						#conceptCol <- colorRampPalette(c('#00ff00',bgColor))                  
						#colorMap <- conceptCol(zeroColorIndex)
						#conceptCol <- colorRampPalette(c(bgColor,'#ff0000'))
						#colorMap <- c(colorMap, conceptCol(colorLevel-zeroColorIndex + 1)) 
						if (symmetry) {
							tempNames <- names(inputValue)
							inputValue <- c(-max(abs(as.numeric(inputValue))), -min(abs(as.numeric(inputValue))), as.numeric(inputValue), 
														min(abs(as.numeric(inputValue))), max(abs(as.numeric(inputValue))))
							names(inputValue) <- c('negativeMax', 'negativeMin', tempNames, 'positiveMin', 'positiveMax')
						}
						zeroColorIndex <- 1 + ceiling(abs(min(as.numeric(inputValue))) * (colorLevel-1) / (max(as.numeric(inputValue))-min(as.numeric(inputValue))))
						conceptCol <- colorRampPalette(c('#3333ff',bgColor))                  
						colorMap <- conceptCol(zeroColorIndex)
						conceptCol <- colorRampPalette(c(bgColor,'#ff3333'))
						colorMap <- c(colorMap, conceptCol(colorLevel-zeroColorIndex + 1))
					}
				} 
			}
			colorValues <- .colorMatch(as.numeric(inputValue), colorMap=colorMap, matchMode=matchMode, zeroColorIndex=zeroColorIndex)
			names(colorValues) <- names(inputValue)
			if (symmetry) {
				inputValue <- rev(rev(inputValue)[-1:-2])[-1:-2]
				colorValues <- rev(rev(colorValues)[-1:-2])[-1:-2]
			}
			V(g)[V(g)$name %in% names(colorValues)]$color <- colorValues[intersect(V(g)$name, names(colorValues))]
   	} else {
			stop('input values do not match input list or input values are not numeric!!!')
		}
	}
	if (is.character(centroidSize) & length(centroidSize) == 1) {
		if (toupper(centroidSize) %in% "GENENUM") V(g)[1:lengthOfRoots]$size <- 100 * degree(g)[1:lengthOfRoots] / sum(degree(g)[1:lengthOfRoots])
		else stop('can not understand!!!')
			
	} else {
		if (is.numeric(centroidSize) & (length(centroidSize) == lengthOfRoots)) {
			V(g)[1:lengthOfRoots]$size <- 180 * (centroidSize[1:lengthOfRoots]) / sum(centroidSize[1:lengthOfRoots])
		} else {
			stop('Size of junction nodes should be numeric or centroidSize should match inputList!')
		}
	}
	V(g)[1:lengthOfRoots]$color <- "#ffd900"
	# apply(get.edges(g13, 1:ecount(g13)), 2, function(x,y) return(y[x]), V(g13)$name)
	#edgesMatrix <- as.data.frame(apply(get.edges(g,1:ecount(g)), 2, function(x,y) return(y[x]), V(g)$name), stringsAsFactors =FALSE)
	#colnames(edgesMatrix) <- c('NODES1', 'NODES2')
	if (output == 'none') {
		return(g)
		#return(c('graph'=list(g), 'vertex.attributes'=list(as.data.frame(cbind('NODES'=V(g)$label, 'NODE_FILL_COLOR'=V(g)$color, 'NODE_SIZE'=3*V(g)$size), stringsAsFactors =FALSE)), 
		#				'edge.attributes'=list(as.data.frame(cbind(edgesMatrix, 'EDGE_COLOR'=E(g)$color, 'EDGE_LINE_WIDTH'=E(g)$width), stringsAsFactors =FALSE))))
	} else {
		if (output == 'fixed') {
			plot.igraph(g, vertex.label.font=2, vertex.label=V(g)$label, vertex.label.color='#666666', vertex.label.cex=1.5, vertex.frame.color=V(g)$color, edge.loop.angle=-90, layout=layout.fruchterman.reingold) 
		} else {
			tkplot(g, vertex.label.font=2, vertex.label=V(g)$label, vertex.label.color='#666666', vertex.label.cex=1.5, vertex.frame.color=V(g)$color, edge.loop.angle=-90, layout=layout.fruchterman.reingold) 
		}
		return(invisible(g))
		#return(invisible(c('graph'=list(g), 'vertex.attributes'=list(as.data.frame(cbind('NODES'=V(g)$label, 'NODE_FILL_COLOR'=V(g)$color, 'NODE_SIZE'=3*V(g)$size), stringsAsFactors =FALSE)), 
		#				'edge.attributes'=list(as.data.frame(cbind(edgesMatrix, 'EDGE_COLOR'=E(g)$color, 'EDGE_LINE_WIDTH'=E(g)$width), stringsAsFactors =FALSE)))))
	}
}

