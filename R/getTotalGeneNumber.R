`getTotalGeneNumber` <- 
function(categoryType=c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'DOLITE', 'KEGG', 'REACTOME.PATH', 'CABIO.PATH'), known=TRUE, annotationLib=c('org.Ag.eg.db', 'org.Bt.eg.db', 'org.Ce.eg.db', 'org.Cf.eg.db', 'org.Dm.eg.db', 
							'org.Dr.eg.db', 'org.EcK12.eg.db', 'org.EcSakai.eg.db', 'org.Gg.eg.db', 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Mmu.eg.db', 'org.Pt.eg.db', 'org.Rn.eg.db', 
							'org.Ss.eg.db', 'org.Xl.eg.db', 'org.At.tair.db', 'org.Pf.plasmo.db', 'org.Sc.sgd.db')) {
	categoryType <- toupper(categoryType)
	annotationLib <- match.arg(annotationLib) 
	require(annotationLib, character.only=TRUE)
	libname <- sub('\\.db', '', annotationLib)
	idType <- switch(sub('org.*[:.:]', '', libname), 'eg'='EG', 'tair'='TAIR', 'ORF')
	if (!known) return(count.mappedkeys(get(paste(libname, 'GENENAME', sep=''))))
	categoryType <- match.arg(categoryType)
	if (categoryType == 'DOLITE') {
	   if (annotationLib == 'org.Hs.eg.db') return(4051)
	   else stop('DOLite is only designed for human! No data for other species ...')
	}
	
	if (categoryType == 'CABIO.PATH') {
		return(switch(annotationLib,
			'org.Hs.eg.db'=.getcaBIOIDInfo('http://cabioapi.nci.nih.gov/cabio43/GetXML?query=gov.nih.nci.cabio.domain.Gene&gov.nih.nci.cabio.domain.Pathway[@name=h_*]', IDinfo='geneNumber'), 
			'org.Mm.eg.db'=.getcaBIOIDInfo('http://cabioapi.nci.nih.gov/cabio43/GetXML?query=gov.nih.nci.cabio.domain.Gene&gov.nih.nci.cabio.domain.Pathway[@name=m_*]', IDinfo='geneNumber'),  
			NULL))
	}
	
	if (categoryType == 'REACTOME.PATH') {
		return(switch(annotationLib,
			'org.At.tair.db'=516, 'org.Ce.eg.db'=627, 'org.Dm.eg.db'=686, 'org.EcK12.eg.db'=185, 'org.EcSakai.eg.db'=185, 'org.Gg.eg.db'=840, 'org.Hs.eg.db'=1019, 'org.Mm.eg.db'=900,  
			'org.Pf.plasmo.db'=308, 'org.Rn.eg.db'=883, 'org.Sc.sgd.db'=473, NULL))
	}
	
	return(switch(categoryType,
		'GO' = count.mappedkeys(get(paste(libname, 'GO', sep=''))),
		'KEGG' = count.mappedkeys(get(paste(libname, 'PATH', sep=''))),
		'GO.BP' = length(unique(unlist(lookUp('GO:0008150', libname, paste("GO2ALL", idType, 'S', sep=''))))),
		'GO.MF' = length(unique(unlist(lookUp('GO:0003674', libname, paste("GO2ALL", idType, 'S', sep=''))))),
		'GO.CC' = length(unique(unlist(lookUp('GO:0005575', libname, paste("GO2ALL", idType, 'S', sep=''))))) ))

}