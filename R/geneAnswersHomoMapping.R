`geneAnswersHomoMapping` <-
function(x, species=c('human', 'rat', 'mouse', 'fly'), speciesL=c('human', 'rat', 'mouse', 'fly'), 
									mappingMethod=c("direct", "biomaRt", "none"), filterGenes=NULL, verbose=TRUE) {
	y <- x
	species <- match.arg(species)
	speciesL <- match.arg(speciesL)
	print('Change annLib ...')
	y@annLib = switch(speciesL,
		'human'='org.Hs.eg.db',
		'rat'='org.Rn.eg.db',
		'mouse'='org.Mm.eg.db',
		'fly'='org.Dm.eg.db')
		
	if (verbose) print('Converting geneInput ...')
	temp <- getHomoGeneIDs(x@geneInput[,1], species=species, speciesL=speciesL, mappingMethod=mappingMethod)
	if (!is.null(filterGenes)) temp <- temp[temp %in% filterGenes]
	y@geneInput <- x@geneInput[x@geneInput[,1] %in% names(temp),]
	y@geneInput[,1] <- temp[y@geneInput[,1]]
	
	if (verbose) print('Converting genesInCategory ...')
	temp <- getHomoGeneIDs(unique(unlist(x@genesInCategory)), species=species, speciesL=speciesL, mappingMethod=mappingMethod) 
	filtConvert <- function(inputGenes, filter, homoIDs) {
		tempGenes <- homoIDs[names(homoIDs) %in% inputGenes]
		if (!is.null(filter)) tempGenes <- intersect(tempGenes, filter)
		if (length(tempGenes) < 1) stop('No homo gene in category! Check convert method')
		return(tempGenes)
	}
	y@genesInCategory <- lapply(x@genesInCategory, filtConvert, filterGenes, temp)
	
	if (!is.null(x@geneExprProfile)) {
		if (verbose) print('Converting geneExprProfile ...')
		temp <- getHomoGeneIDs(rownames(x@geneExprProfile), species=species, speciesL=speciesL, mappingMethod=mappingMethod)
		if (!is.null(filterGenes)) temp <- temp[temp %in% filterGenes]
		y@geneExprProfile <- x@geneExprProfile[rownames(x@geneExprProfile) %in% names(temp),]
		rownames(y@geneExprProfile) <- temp[rownames(y@geneExprProfile)] 
	}
	return(y)
}

