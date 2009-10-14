`geneAnswersBuilder` <-
function(geneInput, annotationLib, categoryType=NULL, testType=c('hyperG', 'none'), totalGeneNumber=NULL, geneExpressionProfile=NULL, categorySubsetIDs=NULL, pvalueT=0.01, FDR.correction=FALSE, 
								verbose=TRUE, sortBy=c('pvalue', 'geneNum', 'foldChange', 'oddsRatio', 'correctedPvalue', 'none'), ...) {
	testType <- match.arg(testType)
	sortBy <- match.arg(sortBy)
	if(!(is.null(totalGeneNumber))) {
		if (!(is.numeric(totalGeneNumber)) & !(tolower(totalGeneNumber) %in% c('human', 'mouse', 'rat', 'fly'))) stop('TotalGeneNumber should be NULL or numeric or one of "human", "mouse", "rat" and "fly"! Abort GeneAnswers Building ...')
	}
	x <- new("GeneAnswers")
	require(annotate)
	if (is.vector(geneInput) | is.data.frame(geneInput) | is.matrix(geneInput)) {
		if (is.vector(geneInput)) {
			geneIDs <- as.character(geneInput)
		} else { 
			geneIDs <- as.character(geneInput[,1])
		}
		x@geneInput <- as.data.frame(geneInput)
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
		if (is.character(annotationLib) & (length(annotationLib) == 1) & (length(categoryType) == 1)) {
			if (toupper(categoryType) %in% c('GO','GO.BP','GO.CC','GO.MF','DOLITE','KEGG')) {
				if (length(which(c('org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db', 'org.Dm.eg.db') %in% annotationLib)) == 0) {
					stop(paste(annotationLib, 'can not be loaded! Abort GeneAnswers Building ...'))
				} else {
					switch(which(c('org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db', 'org.Dm.eg.db') %in% annotationLib),
						require('org.Hs.eg.db'),
						require('org.Mm.eg.db'),
						require('org.Rn.eg.db'),
						require('org.Dm.eg.db'))
				}
				data('DOLite', package='GeneAnswers')
   				annLibList = switch(toupper(categoryType), 
					'GO'=getGOList(geneIDs, annotationLib, GOCat='ALL', ...),
					'GO.BP'=getGOList(geneIDs, annotationLib, GOCat='BP', ...),
					'GO.CC'=getGOList(geneIDs, annotationLib, GOCat='CC', ...),
					'GO.MF'=getGOList(geneIDs, annotationLib, GOCat='MF', ...),
					'DOLITE'=DOLite,
					'KEGG'= getPATHList(geneIDs, annotationLib))
					x@annLib <- annotationLib
					x@categoryType <- categoryType
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

	x@genesInCategory <- list()
	for (i in 1:length(testLibList)) {
		temp <- list(intersect(testLibList[[i]], geneIDs))
		if (length(temp[[1]]) > 0) {
			names(temp) <- names(testLibList[i])
			x@genesInCategory <- c(x@genesInCategory, temp)
		}
	}
	if (verbose) print('genesInCategory has built in ...')

	if (testType == 'hyperG') {
		if (is.null(totalGeneNumber)) {
			if(is.null(x@annLib)) stop('Missing total gene number for hypergeometic test! Abort GeneAnswers Building ...')
			if (x@annLib %in% c('org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db', 'org.Dm.eg.db')) {
				totalGeneNumber <- switch(x@annLib,
					'org.Hs.eg.db'=45384,
					'org.Mm.eg.db'=61498,
					'org.Rn.eg.db'=37536,
					'org.Dm.eg.db'=22606)
			} else stop('Missing total gene number for hypergeometic test! Abort GeneAnswers Building ...')
		} else {
			if (tolower(totalGeneNumber) %in% c('human', 'mouse', 'rat', 'fly')) {
				totalGeneNumber <- switch(totalGeneNumber,
					'human'=45384,
					'mouse'=61498,
					'rat'=37536,
					'fly'=22606)
			}
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
		fullResult <- fullResult[fullResult[, 'fdr p value'] < pvalueT,]
	} else {
		fullResult <- fullResult[fullResult[, 'p value'] < pvalueT, 1:(dim(fullResult)[2]-1)]
	}
	x@enrichmentInfo <- fullResult
#	x@genesInCategory <- x@genesInCategory[rownames(x@enrichmentInfo)] 

	if (verbose) print('testType, pvalueT and enrichmentInfo have built in ...')
	if (is.null(geneExpressionProfile)) x@geneExprProfile <- NULL
	else x@geneExprProfile <- as.data.frame(geneExpressionProfile)
	if (verbose) print('geneExpressionProfile has been built in ...')
	
	print('GeneAnswers instance has been successfully generated!')
	if (sortBy != 'none') return(geneAnswersSort(x, sortBy=sortBy))
	else return(x)
}

