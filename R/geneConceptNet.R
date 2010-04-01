`geneConceptNet` <-
function(inputList, lengthOfRoots=NULL, inputValue=NULL, centroidSize='geneNum', output=c('fixed','interactive'), colorMap=NULL, matchMode=c('absolute', 'relative'), zeroColorIndex=NULL) {
	require(igraph)
	output <- match.arg(output)
	matchMode <- match.arg(matchMode)
	
	if (is.null(lengthOfRoots)) lengthOfRoots <- length(inputList)
	temp <- .list2graph(inputList)
	g <- temp[[1]]
	nodes <- temp[[2]]
	vertexLabels <- names(nodes)
	vertexLabels[1:length(inputList)] <- names(inputList) 
	E(g)$color <- "#9999ff"
	E(g)$width <- 3
	V(g)$size <- 8
	V(g)$color <- "#dddddd"
	V(g)$label <- vertexLabels
	if (!is.null(inputValue)) {
		if(all(!is.na(as.numeric(inputValue))) & (all(unique(unlist(inputList)) %in% names(inputValue)))) {
			if (is.null(colorMap)) {
				colorLevel <- 256
				if ((min(as.numeric(inputValue)) > 0) | (max(as.numeric(inputValue)) < 0)) {
					if (min(as.numeric(inputValue)) > 0) {
						zeroColorIndex <- 1
						conceptCol <- colorRampPalette(c('#ffffff','#ff0000'))
					} else {
						zeroColorIndex <- colorLevel
						conceptCol <- colorRampPalette(c('#00ff00','#ffffff'))
					}
					colorMap <- conceptCol(colorLevel)
				} else {
					zeroColorIndex <- 1 + ceiling(abs(min(as.numeric(inputValue))) * (colorLevel-1) / (max(as.numeric(inputValue))-min(as.numeric(inputValue))))
					conceptCol <- colorRampPalette(c('#00ff00','#ffffff'))                  
					colorMap <- conceptCol(zeroColorIndex)
					conceptCol <- colorRampPalette(c('#ffffff','#ff0000'))
					colorMap <- c(colorMap, conceptCol(colorLevel-zeroColorIndex + 1))
				}
			}
			colorValues <- .colorMatch(as.numeric(inputValue), colorMap=colorMap, matchMode=matchMode, zeroColorIndex=zeroColorIndex)
			names(colorValues) <- names(inputValue)
			V(g)[nodes[intersect(names(nodes), names(colorValues))]]$color <- colorValues[intersect(names(nodes), names(colorValues))]
   	} else {
			stop('input values do not match input list or input values are not numeric!!!')
		}
	}
	if (is.character(centroidSize) & length(centroidSize) == 1) {
		if (toupper(centroidSize) %in% "GENENUM") V(g)[0:(lengthOfRoots-1)]$size <- 100 * degree(g)[1:lengthOfRoots] / sum(degree(g)[1:lengthOfRoots])
		else stop('can not understand!!!')
			
	} else {
		if (is.numeric(centroidSize) & (length(centroidSize) == lengthOfRoots)) {
			V(g)[0:(lengthOfRoots-1)]$size <- 180 * (centroidSize[1:lengthOfRoots]) / sum(centroidSize[1:lengthOfRoots])
		} else {
			stop('Size of junction nodes should be numeric or centroidSize should match inputList!')
		}
	}
	V(g)[0:(lengthOfRoots-1)]$color <- "#ffd900"
	if (output == 'fixed') {
		plot.igraph(g, vertex.label.font=2, vertex.label=V(g)$label, vertex.label.color='#666666', vertex.label.cex=1.5, vertex.frame.color=V(g)$color, edge.loop.angle=-90, layout=layout.fruchterman.reingold) 
	} else {
		tkplot(g, vertex.label.font=2, vertex.label=V(g)$label, vertex.label.color='#666666', vertex.label.cex=1.5, vertex.frame.color=V(g)$color, edge.loop.angle=-90, layout=layout.fruchterman.reingold) 
	}
	
}

