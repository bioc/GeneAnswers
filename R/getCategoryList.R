`getCategoryList` <-
function(geneVector, lib, categoryType) {
	if (lib %in% c('org.Bt.eg.db', 'org.Ce.eg.db', 'org.Cf.eg.edu', 'org.Dm.eg.db', 
					'org.Dr.eg.db', 'org.EcK12.eg.db', 'org.EcSakai.eg.db', 'org.Gg.eg.db', 
					'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db', 'org.Ss.eg.db')) {
						require(lib, character.only=TRUE)
						libname <- paste(sub('.db', '', lib), categoryType, sep='')
						categoryTerms <- unique(unlist(lookUp(geneVector, lib, categoryType)))
						if (!is.null(categoryTerms)){
							return(lookUp(categoryTerms[!(categoryTerms %in% NA)], lib, paste(categoryType, '2EG', sep='')))
						}
						else stop('input genes do not belong to any category!!!') 
					} else {
						stop('the specified annotation library can not be found!!!')
					}
}

