`geneConceptNet` <-
function(inputList, inputValue=NULL, centroidSize='geneNum', output=c('fixed','interactive')) {
	require(igraph)
	output <- match.arg(output)
	
	temp <- .list2graph(inputList, saveFile=FALSE)
	nodes <- temp[[2]]
	vertexLabels <- names(nodes)
	vertexLabels[1:length(inputList)] <- names(inputList) 
	g <- graph.edgelist(temp[[1]], directed=FALSE)
	E(g)$color <- "#8888ff"
	V(g)$size <- 8
	V(g)$color <- "#dddddd"
	V(g)$label <- vertexLabels
	if (!is.null(inputValue)) {
		if(is.numeric(inputValue) & (all(unique(unlist(inputList)) %in% names(inputValue)))) {
			V(g)[nodes[names(nodes) %in% names(inputValue[inputValue < 0])]]$color <- "#00ff00"
			V(g)[nodes[names(nodes) %in% names(inputValue[inputValue > 0])]]$color <- "#ff0000"
		} else {
			stop('input value does not match input list!!!')
		}
	}
	if (is.character(centroidSize) & length(centroidSize) == 1) {
		if (toupper(centroidSize) %in% "GENENUM") V(g)[0:(length(inputList)-1)]$size <- 100 * degree(g)[1:length(inputList)] / sum(degree(g)[1:length(inputList)])
		else stop('can not understand!!!')
			
	} else {
		if (is.numeric(centroidSize) & (length(centroidSize) == length(inputList))) {
			V(g)[0:(length(inputList)-1)]$size <- 180 * (centroidSize[1:length(inputList)]) / sum(centroidSize[1:length(inputList)])
		} else {
			stop('Size of junction nodes should be numeric or centroidSize should match inputList!')
		}
	}
	V(g)[0:(length(inputList)-1)]$color <- "#ffd900"
	if (output == 'fixed') {
		plot(g, vertex.label.font=2, vertex.label=V(g)$label, vertex.label.color='#666666', vertex.label.cex=1.5, layout=layout.fruchterman.reingold) 
	} else {
		tkplot(g, vertex.label.font=2, vertex.label=V(g)$label, vertex.label.color='#666666', vertex.label.cex=1.5, layout=layout.fruchterman.reingold) 
	}
	
}

