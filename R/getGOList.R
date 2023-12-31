`getGOList` <-
function(geneVector, lib, GOCat = c('ALL', 'BP', 'CC', 'MF'), level = 1) {
	GOCat <- match.arg(GOCat)
	if (lib %in% c( 'org.Ag.eg.db', 'org.Bt.eg.db', 'org.Ce.eg.db', 'org.Cf.eg.db', 'org.Dm.eg.db', 'org.Dr.eg.db', 'org.EcK12.eg.db', 'org.EcSakai.eg.db', 'org.Gg.eg.db', 
					'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Mmu.eg.db', 'org.Pt.eg.db', 'org.Rn.eg.db', 'org.Ss.eg.db', 'org.Xl.eg.db', 
					'org.At.tair.db', 'org.Pf.plasmo.db', 'org.Sc.sgd.db')) {
						require(lib, character.only=TRUE)
						require(GO.db)
						idType <- switch(sub('org.*[:.:]', '', sub('\\.db', '', lib)),
									'eg'='EG', 'tair'='TAIR', 'ORF')
						goIDs <- unique(unlist(lapply(lookUp(geneVector, lib, 'GO'), names)))
						if (GOCat != 'ALL' ) goIDs <- goIDs[goIDs %in% aqListGOIDs(GOCat)]
						allGOIDs <- goIDs
						if ((GOCat == 'ALL') | (GOCat == 'BP'))	allGOIDs <- c(allGOIDs, unique(unlist(lookUp(goIDs, "GO", "BPANCESTOR")))) 
						if ((GOCat == 'ALL') | (GOCat == 'CC'))	allGOIDs <- c(allGOIDs, unique(unlist(lookUp(goIDs, "GO", "CCANCESTOR"))))
						if ((GOCat == 'ALL') | (GOCat == 'MF')) allGOIDs <- c(allGOIDs, unique(unlist(lookUp(goIDs, "GO", "MFANCESTOR"))))
						allGOIDs <- unique(allGOIDs[!(allGOIDs %in% NA)])
						naGOIDs <- unlist(lapply(lookUp(allGOIDs, lib, paste("GO2ALL", idType, 'S', sep='')), is.na))
						allGOIDs <- allGOIDs[!(allGOIDs %in% names(naGOIDs[naGOIDs %in% TRUE]))]
						allGOIDs <- allGOIDs[!(allGOIDs %in% .filterGOIDs(GOCategory=GOCat, level=level))]
						if (!is.null(allGOIDs)) {
							return(lapply(lookUp(allGOIDs, lib, paste("GO2ALL", idType, 'S', sep='')), unique))
						} else stop('input genes do not belong to any category!!!') 
	} else {
		stop('the specified annotation library can not be found!!!')
	}
}

