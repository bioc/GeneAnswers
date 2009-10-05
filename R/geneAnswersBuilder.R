`geneAnswersBuilder` <-
function(geneInput, annotationLib, categoryType=NULL, testType=c('hyperG', 'none'), geneExpressionProfile=NULL, 
								categorySubsetIDs=NULL, pvalueT=0.01, FDR.correction=FALSE, verbose=TRUE, ...) {
	testType <- match.arg(testType)
	x <- new("GeneAnswers")
	require(annotate)
	if (is.vector(geneInput) | is.data.frame(geneInput) | is.matrix(geneInput)) {
		if (is.vector(geneInput)) {
			geneIDs <- geneInput
		} else { 
				geneIDs <- geneInput[,1]
		}
		x@geneInput <- as.data.frame(geneInput)
	} else {
		stop('input gene info can not be recognized!')
	}
	if (verbose) print('geneInput has built in ...')

	if (is.list(annotationLib) & (is.null(categoryType))){
		annLibList <- annotationLib
		x@annLib <- NULL
		x@categoryType <- 'User defiend'
	} else {
		if (is.character(annotationLib) & (length(annotationLib) == 1)) {
			if ((length(categoryType) == 1)) {
				annLibList = switch(toupper(categoryType), 
					'GO'=getGOList(geneIDs, annotationLib, GOCat='ALL', ...),
					'GO.BP'=getGOList(geneIDs, annotationLib, GOCat='BP', ...),
					'GO.CC'=getGOList(geneIDs, annotationLib, GOCat='CC', ...),
					'GO.MF'=getGOList(geneIDs, annotationLib, GOCat='MF', ...),
					'DOLITE'=DOLite,
					'KEGG'= getPATHList(geneIDs, annotationLib))
					x@annLib <- annotationLib
					x@categoryType <- categoryType
			}
		} else {
			stop('annotation library can not be recognized!!!')
		}
	}
	if (!is.null(categorySubsetIDs)) {
		testLibList <- annLibList[names(annLibList) %in% categorySubsetIDs]
		if (is.null(testLibList)) stop('test category list is empty!!!')
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

	if (testType == 'hyperG') fullResult <- .hyperGTest(geneIDs, testLibList)
	x@testType <- testType
	if (is.numeric(pvalueT)) {
		if (pvalueT <= 0 | pvalueT > 1) {
			stop('pvalue threshold should be a number in (0,1]!')
		}
	} else {
		stop('pvalue threshold should be a number in (0,1]!')
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
	else x@geneExprProfile <- geneExpressionProfile
	if (verbose) print('geneExpressionProfile has been built in ...')
	
	print('GeneAnswers instance has been successfully generated!')
	return(x)
}

