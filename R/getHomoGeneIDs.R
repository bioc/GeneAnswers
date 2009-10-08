`getHomoGeneIDs` <-
function (oriGeneIDs, species=c('human', 'rat', 'mouse', 'yeast', 'fly'), speciesL=c('human', 'rat', 'mouse', 'yeast', 'fly'), 
							mappingMethod=c("direct", "biomaRt", "none")) {
	species <- match.arg(species)
	speciesL <- match.arg(speciesL)
	mappingMethod <- match.arg(mappingMethod)
	if (NA %in% as.numeric(oriGeneIDs)) {
		print('Only EntrezIDs are supported and non-EntrezIDs will be removed!')
		temp <- as.numeric(oriGeneIDs)
		temp <- as.character(temp[!(temp %in% NA)])
		oriGeneIDs <- oriGeneIDs[oriGeneIDs %in% temp]
		}
	if (mappingMethod == 'direct') {
		#print('This method currently only work between human and mouse!')
		if ((species != 'human') & (species != 'mouse') & (speciesL != 'human') & (speciesL != 'mouse')) stop('Direct method does not support specified species!')
		lib = switch(species,
			'human'='org.Hs.eg.db',
 			'mouse'='org.Mm.eg.db')
		libL = switch(speciesL,
			'human'='org.Hs.eg.db',
 			'mouse'='org.Mm.eg.db')	
		require(lib, character.only=TRUE)
		require(libL, character.only=TRUE)
		geneSymbols <-unlist(getSymbols(oriGeneIDs, lib, missing='name'))
		if (speciesL == 'human') newGeneSymbols <- toupper(geneSymbols)
		else newGeneSymbols <- paste(toupper(substring(geneSymbols,1,1)), {tolower(substring(geneSymbols,2))}, sep = "")
		names(newGeneSymbols) <- names(geneSymbols)
		newEntrezIDs <- unlist(lookUp(newGeneSymbols, libL, 'SYMBOL2EG'))
		newEntrezIDs <- newEntrezIDs[newGeneSymbols]
		tempV <- names(newGeneSymbols)
		names(tempV) <- newGeneSymbols
		names(newEntrezIDs) <- tempV[names(newEntrezIDs)]
		homogenes <- newEntrezIDs[!(newEntrezIDs %in% NA)]
		if (length(oriGeneIDs) > length(unique(names(homogenes)))) print('Warning: homogenes of some input genes can not be found and are removed!!!') 
		return(homogenes)
	}
	if (mappingMethod == 'biomaRt') {
		require(biomaRt)
		lib = switch(species,
			'human'='hsapiens_gene_ensembl',
	  		'rat'='rnorvegicus_gene_ensembl',
	  		'mouse'='mmusculus_gene_ensembl',
	  		'yeast'='scerevisiae_gene_ensembl',
	  		'fly'='dmelanogaster_gene_ensembl')
		libL = switch(speciesL,
			'human'='hsapiens_gene_ensembl',
	  		'rat'='rnorvegicus_gene_ensembl',
	  		'mouse'='mmusculus_gene_ensembl',
	  		'yeast'='scerevisiae_gene_ensembl',
	  		'fly'='dmelanogaster_gene_ensembl')
		from = useMart("ensembl", dataset = lib)
		to = useMart("ensembl", dataset = libL)
		print('Converting ...') 
		homogenesM <- as.matrix(getLDS(attributes = c("entrezgene"), filters = "entrezgene", values = oriGeneIDs, mart = from,  attributesL = c("entrezgene"), martL = to))
		homogenes <- as.character(homogenesM[,2])
		names(homogenes) <- homogenesM[,1]
		if (length(oriGeneIDs) > length(unique(names(homogenes)))) print('Warning: homogenes of some input genes can not be found and are removed!!!')
		sortHomo <- function(input, valuesVector) { temp <- valuesVector[names(valuesVector) %in% input]; return(temp[input[input %in% names(valuesVector)]])}
		return(unlist(lapply(oriGeneIDs[oriGeneIDs %in% unique(names(homogenes))], sortHomo, homogenes)))
	}
}

