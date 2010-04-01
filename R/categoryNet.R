`categoryNet` <- 
function(catGenesList, centroidSize=NULL, output=c('fixed','interactive')) {
	output <- match.arg(output)  
	nodes <- c(0:(length(catGenesList)-1))
	names(nodes) <- names(catGenesList)
	edgeMatrix <- c()
	edgeWidth <- c()
	monoNodes <- c(1:length(nodes)) - 1
	for (i in 1:(length(catGenesList) - 1)) {
		for (j in (i+1):length(catGenesList)) {
			if (length(intersect(catGenesList[[i]], catGenesList[[j]])) > 0) {
				edgeMatrix <- c(edgeMatrix, c((i-1), (j-1)))
				edgeWidth <- c(edgeWidth, length(intersect(catGenesList[[i]], catGenesList[[j]])))
			} 
		}
	}
	
	g <- igraph::graph(edgeMatrix, n = length(nodes), direct = FALSE)
	E(g)$color <- "#9999ff"
	widthValues <- edgeWidth/length(unique(unlist(catGenesList)))
	scaledWidth <- scale(widthValues, scale=(max(widthValues)-min(widthValues)))
	E(g)$width <- as.integer((scaledWidth - min(scaledWidth)) * 11) + 1
	E(g)$label <- as.character(edgeWidth)
	V(g)$color <- "#ffd900"
	V(g)$label <- names(nodes)
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