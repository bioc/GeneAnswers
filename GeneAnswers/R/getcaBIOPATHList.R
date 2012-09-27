`getcaBIOPATHList` <- 
function(lls) {
	genes <- unique(unlist(entrez2caBIO(unlist(lls))))
	if (is.null(genes)) stop('No caBIO genes mapping to the given Entrez genes! Aborting ...')
	pathways <- unique(unlist(.getcaBIOPATHs(genes)))
	if (is.null(pathways)) stop('No pathways in caBIO found! Aborting ...')
	return(.mapListcaBIOIDs2entrez(.getcaBIOPATHGenes(pathways)))
}