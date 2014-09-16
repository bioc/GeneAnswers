`geneAnswersBuilder` <-
function(geneInput, annotationLib, categoryType=NULL, testType=c('hyperG', 'none'), known=TRUE, totalGeneNumber=NULL, geneExpressionProfile=NULL, categorySubsetIDs=NULL, pvalueT=0.01, FDR.correction=FALSE, 
								verbose=TRUE, sortBy=c('pvalue', 'geneNum', 'foldChange', 'oddsRatio', 'correctedPvalue', 'none'), ...) {
	testType <- match.arg(testType)
	sortBy <- match.arg(sortBy)
	if(!(is.null(totalGeneNumber))) {
		if (!(is.numeric(totalGeneNumber)) & !(tolower(totalGeneNumber) %in% c('anopheles', 'arabidopsis', 'bovine', 'worm', 'canine', 'fly', 'zebrafish', 'ecolistraink12', 'ecolistrainsakai', 'chicken', 
																				'human', 'mouse', 'rhesus', 'malaria', 'chimp', 'rat', 'yeast', 'pig', 'xenopus') ))
			stop('TotalGeneNumber should be NULL or numeric or one of anopheles, arabidopsis, bovine, worm, canine, fly, zebrafish, ecolistraink12, ecolistrainsakai, chicken, human, mouse, rhesus, malaria, chimp, rat, yeast, pig and xenopus! Abort GeneAnswers Building ...')
	}
	x <- new("GeneAnswers")
	if (is.vector(geneInput) | is.data.frame(geneInput) | is.matrix(geneInput)) {
		if (is.vector(geneInput)) {
			geneIDs <- as.character(geneInput)
			geneInput <- geneIDs
		} else { 
			geneIDs <- as.character(geneInput[,1])
			geneInput[,1] <- geneIDs
		}
		if (NA %in% geneIDs) stop('Given gene IDs contain NA, please check!!! Aborting ...')
		x@geneInput <- as.data.frame(geneInput, stringsAsFactors=FALSE)
		#rownames(x@geneInput) <- geneIDs
	} else {
		stop('input gene info can not be recognized! Abort GeneAnswers Building ...')
	}
	if (verbose) print('geneInput has built in ...')

	if (is.list(annotationLib)){
		annLibList <- annotationLib
		x@annLib <- NULL
		if (!(is.null(categoryType))) print('categoryType is set User defined')
		x@categoryType <- 'User defined'
	} else {
		if (is.null(annotationLib) & !(is.null(totalGeneNumber))) {
			annotationLib <- .name2lib(totalGeneNumber)                                                                                                                                                                                      
		}
		if (is.character(annotationLib) & (length(annotationLib) == 1) & (length(categoryType) == 1)) {
			if (toupper(categoryType) %in% c('GO','GO.BP','GO.CC','GO.MF','DOLITE','KEGG', 'REACTOME.PATH', 'CABIO.PATH')) {
				require(annotationLib, character.only=TRUE)
				data('DOLite', package='GeneAnswers')
   				annLibList = switch(toupper(categoryType), 
					'GO'=getGOList(geneIDs, annotationLib, GOCat='ALL', ...),
					'GO.BP'=getGOList(geneIDs, annotationLib, GOCat='BP', ...),
					'GO.CC'=getGOList(geneIDs, annotationLib, GOCat='CC', ...),
					'GO.MF'=getGOList(geneIDs, annotationLib, GOCat='MF', ...),
					'DOLITE'=DOLite,
					'KEGG'= getPATHList(geneIDs, annotationLib),
					'REACTOME.PATH'=getREACTOMEPATHList(geneIDs),
					'CABIO.PATH'=stop("Due to termination of caBig, this function is removed in this version!"))
					x@annLib <- annotationLib
					x@categoryType <- toupper(categoryType)
			} else {
				stop('AnnotationLib can not be recognized! Abort GeneAnswers Building ...')
			}
		} else {
			stop('Annotation library can not be recognized!  Abort GeneAnswers Building ...')
		}
	}
	if (!is.null(categorySubsetIDs)) {
		testLibList <- annLibList[names(annLibList) %in% categorySubsetIDs]
		if (is.null(testLibList)) stop('Test category list is empty! Abort GeneAnswers Building ...')
	} else {
		testLibList <- annLibList
	}
	if (verbose) print('annLib and categoryType have built in ...')

#	x@genesInCategory <- list()
#	for (i in 1:length(testLibList)) {
#		temp <- list(intersect(testLibList[[i]], geneIDs))
#		if (length(temp[[1]]) > 0) {
#			names(temp) <- names(testLibList[i])
#			x@genesInCategory <- c(x@genesInCategory, temp)
#		}
#	}
	
	x@genesInCategory <- lapply(testLibList, function(x,y) return(unique(intersect(x,y))), geneIDs)
	x@genesInCategory <- x@genesInCategory[which(sapply(x@genesInCategory, length) != 0)]
	testLibList <- testLibList[names(x@genesInCategory)]
	if (verbose) print('genesInCategory has built in ...')

	if (testType == 'hyperG') {
		if (known) {
			if (verbose) print('Enrichment test is only performed based on annotated genes')
			geneIDs <- geneIDs[geneIDs %in% unique(unlist(testLibList))]
		}
		if (x@categoryType == 'User defined') {
			if (!is.numeric(totalGeneNumber)) {
				if (known) {
					totalGeneNumber <- length(unique(unlist(testLibList)))
				} else {
					if (is.null(totalGeneNumber)) stop('totalGeneNumber is not specified! Aborting GeneAnswers Building ...')
					annLibName <- .name2lib(species=totalGeneNumber)
					if (is.null(annLibName)) stop('species is not support! Aborting GeneAnswers Building ...')
					totalGeneNumber <- getTotalGeneNumber(categoryType='GO', known=FALSE, annotationLib=annLibName)
				}
			}
		} else {
			if (!is.numeric(totalGeneNumber)) totalGeneNumber <- max(getTotalGeneNumber(categoryType=x@categoryType, known=known, annotationLib=x@annLib), length(unique(unlist(testLibList)))) 
			if (is.null(totalGeneNumber)) stop(paste(x@categoryType,' does not support your specified species currently! Aborting GeneAnswers Building ...', sep=''))
		}
		fullResult <- .hyperGTest(geneIDs, testLibList, totalNGenes=totalGeneNumber) 
	}
	
	x@testType <- testType
	if (is.numeric(pvalueT)) {
		if (pvalueT <= 0 | pvalueT > 1) {
			stop('pvalue threshold should be a number in (0,1]! Abort GeneAnswers Building ...')
		}
	} else {
		stop('pvalue threshold should be a number in (0,1]! Abort GeneAnswers Building ...')
	}
	x@pvalueT <- pvalueT
	if (FDR.correction) {
		fullResult <- fullResult[fullResult[, 'fdr p value'] <= pvalueT,]
	} else {
		fullResult <- fullResult[fullResult[, 'p value'] <= pvalueT, 1:(dim(fullResult)[2]-1)]
	}
	x@enrichmentInfo <- fullResult
#	x@genesInCategory <- x@genesInCategory[rownames(x@enrichmentInfo)] 

	if (verbose) print('testType, pvalueT and enrichmentInfo have built in ...')
	if (is.null(geneExpressionProfile)) x@geneExprProfile <- NULL
	else x@geneExprProfile <- as.data.frame(geneExpressionProfile, stringsAsFactors=FALSE)
	if (verbose) print('geneExpressionProfile has been built in ...')
	
	print('GeneAnswers instance has been successfully created!')
	if (sortBy != 'none') return(geneAnswersSort(x, sortBy=sortBy))
	else return(x)
}

