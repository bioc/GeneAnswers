`getCategoryList` <-
function(geneVector, lib, categoryType) {
	if (lib %in% c('org.Ag.eg.db', 'org.Bt.eg.db', 'org.Ce.eg.db', 'org.Cf.eg.db', 'org.Dm.eg.db', 
				   'org.Dr.eg.db', 'org.EcK12.eg.db', 'org.EcSakai.eg.db', 'org.Gg.eg.db', 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Mmu.eg.db', 'org.Pt.eg.db', 'org.Rn.eg.db', 
				   'org.Ss.eg.db', 'org.Xl.eg.db', 'org.At.tair.db', 'org.Pf.plasmo.db', 'org.Sc.sgd.db')) {
						require(lib, character.only=TRUE)
						libname <- sub('.db', '', lib)
						idType <- switch(sub('org.*[:.:]', '', libname), 'eg'='EG', 'tair'='TAIR', 'ORF') 
						categoryTerms <- unique(unlist(lookUp(geneVector, lib, categoryType)))
						if (!is.null(categoryTerms)){
							return(lookUp(categoryTerms[!(categoryTerms %in% NA)], lib, paste(categoryType, paste('2', idType, sep=''), sep='')))
						}
						else stop('input genes do not belong to any category!!!') 
					} else {
						stop('the specified annotation library can not be found!!!')
					}
}

