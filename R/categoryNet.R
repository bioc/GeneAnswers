`categoryNet` <- 
function(catGenesList, centroidSize=NULL, output=c('fixed','interactive')) {
	output <- match.arg(output)  
	nodes <- c(0:(length(catGenesList)-1))
	names(nodes) <- names(catGenesList)
	edgeMatrix <- c()
	edgeWidth <- c()
	for (i in 1:(length(catGenesList) - 1)) {
		for (j in (i+1):length(catGenesList)) {
			if (length(intersect(catGenesList[[i]], catGenesList[[j]])) > 0) {
				edgeMatrix <- c(edgeMatrix, c((i-1), (j-1)))
				edgeWidth <- c(edgeWidth, length(intersect(catGenesList[[i]], catGenesList[[j]])))
			}
		}
	}
	g <- graph.edgelist(t(matrix(edgeMatrix, nrow=2)), directed=FALSE)
	E(g)$color <- "#6666ff"
	E(g)$width <- as.integer(edgeWidth/sum(edgeWidth) * 100) + 1
	V(g)$color <- "#ffd900"
	V(g)$label <- names(nodes)
	if (is.null(centroidSize)) V(g)$size <- 10
	else {
		if (is.numeric(centroidSize) & (length(centroidSize) == length(catGenesList))) V(g)$size <- centroidSize
		else  stop('Size of junction nodes should be numeric or centroidSize should match inputList!')
	}
	if (output == 'fixed') {
		plot(g, vertex.label.font=2, vertex.label=V(g)$label, vertex.label.color='#666666', vertex.label.cex=1.5, layout=layout.circle) 
	} else {
		tkplot(g, vertex.label.font=2, vertex.label=V(g)$label, vertex.label.color='#666666', vertex.label.cex=1.5, layout=layout.circle) 
	}
}