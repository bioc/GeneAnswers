`categoryNet` <- 
function(catGenesList, centroidSize=NULL, output=c('fixed','interactive')) {
	output <- match.arg(output)  
	if ((length(catGenesList) == 1) | is.null(unique(unlist(catGenesList))) | all(is.na(unique(unlist(catGenesList))))) stop('the given category(s) can not be build a network, Aborting ... !')
	edgesMatrix <- c()
	edgeWidth <- c()
	for (i in 1:(length(catGenesList)-1)) {
		edgesMatrix <- rbind(edgesMatrix, cbind(rep(names(catGenesList)[i], length(names(catGenesList)[-1:-i])), names(catGenesList)[-1:-i])) 
		edgeWidth <- c(edgeWidth, sapply(lapply(catGenesList[c(-1:-i)], intersect, catGenesList[[i]]), length))
	}
	
	# remove the edges and nodes without sharing anything
	connection <- which(edgeWidth > 0)
	edgesMatrix <- edgesMatrix[connection,]
	edgeWidth <- edgeWidth[connection]
	
	g <- igraph::graph.edgelist(edgesMatrix, direct = FALSE)
	connectedNodes <- unique(as.vector(edgesMatrix))
	if (length(connectedNodes) < length(catGenesList)) {
		g <- g + vertices(names(catGenesList)[!(names(catGenesList) %in% connectedNodes)])
	}
	E(g)$color <- "#9999ff"
	widthValues <- edgeWidth/length(unique(unlist(connectedNodes)))
	scaledWidth <- scale(widthValues, scale=(max(widthValues)-min(widthValues)))
	E(g)$width <- as.integer((scaledWidth - min(scaledWidth)) * 11) + 1
	E(g)$label <- as.character(edgeWidth)
	V(g)$color <- "#ffd900"
	V(g)$label <- V(g)$name
	if (is.null(centroidSize)) V(g)$size <- 10
	else {
		if (is.numeric(centroidSize) & (length(centroidSize) == length(catGenesList))) V(g)$size <- centroidSize
		else  stop('Size of junction nodes should be numeric or centroidSize should match inputList!')
	}
	if (output == 'fixed') {
		plot.igraph(g, vertex.label.font=2, vertex.label=V(g)$label, vertex.label.color='#666666', vertex.label.cex=1.5, edge.label.font=4, 
				edge.label=E(g)$label, edge.label.color='#ff5500', edge.label.cex=1.5, layout=layout.circle) 
	} else {
		tkplot(g, vertex.label.font=2, vertex.label=V(g)$label, vertex.label.color='#666666', vertex.label.cex=1.5, edge.label.font=4, 
				edge.label=E(g)$label, edge.label.color='#ff5500', edge.label.cex=1.5, layout=layout.circle) 
	}
}