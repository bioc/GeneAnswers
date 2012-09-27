geneFunSummarize <- function(genes, gene2Onto, Onto2offspring, rmOntoID=c("DOID:4", "DOID:63"), p.value.th=0.01, fdr.adjust='fdr', minNumTh=2, includeTestOnto=TRUE, directMapConstraint=FALSE) {
	
	numAllEvd <- length(unlist(gene2Onto))
	
	## convert the graph as a list of edges
	if (class(Onto2offspring) == 'graphNEL') {
		Onto2offspring <- edges(Onto2offspring)
	}
	OntoIDwithOffspring <- names(Onto2offspring)

	## Produce a reverse map from Ontology to gene (direct evidence based on the mapping from gene 2 ontology
	mapLen <- sapply(gene2Onto, length)
	geneNonto <- data.frame(gene=rep(names(gene2Onto), mapLen), Ontology=unlist(gene2Onto))
	Onto2directEvd <- tapply(geneNonto[, 1], geneNonto[, 2], function(x) list(x))
	num.Onto2directEvd <- sapply(Onto2directEvd, length)
	
	## Estimate the total number of evidences of each ontology including its offspring  
	num.Onto2allEvd <- sapply(Onto2offspring, function(x) {
		x <- x[x %in% names(num.Onto2directEvd)]
		sum(num.Onto2directEvd[x])
	})

	## Include the direct evidence of the terms themselves
	if (includeTestOnto) {
		num.Onto2allEvd[names(num.Onto2directEvd)] <- num.Onto2allEvd[names(num.Onto2directEvd)] + num.Onto2directEvd
	}
	
	if (is.null(genes)) genes <- names(gene2Onto)
	allDirectMapID <- names(Onto2directEvd)
	
	testSigInfo <- NULL
	for (gene.i in genes) {
		testSigInfo.i <- NULL		# keep the information of how the test was performed (the evidence and p-values of related ontologies)
		allSigOntoID.i <- NULL
		
		relatedOntoID.direct.i <- gene2Onto[[gene.i]]
		if (length(relatedOntoID.direct.i) < minNumTh) {
			testSigInfo <- c(testSigInfo, list(list(allEvidence=relatedOntoID.direct.i)))		
			next
		}
		relatedOntoID.i <- OntoIDwithOffspring[sapply(Onto2offspring, function(x) sum(relatedOntoID.direct.i %in% x) >= minNumTh)]
		if (directMapConstraint) {
			relatedOntoID.i <- intersect(relatedOntoID.i, allDirectMapID)
		}
		relatedOntoID.i <- unique(c(relatedOntoID.i, relatedOntoID.direct.i))

		# remove the predefined Ontology ids (usually becasue they are too general)
		if (!is.null(rmOntoID))	relatedOntoID.i <- relatedOntoID.i[!(relatedOntoID.i %in% rmOntoID)]
		if (length(relatedOntoID.i) == 0) {
			testSigInfo <- c(testSigInfo, list(list(allEvidence=relatedOntoID.direct.i)))		
			next
		}

		# The test will go through all relatedOntoID.i
		bestOntoID.i <- bestPvalue.i <- bestP.adjust.i <- NULL

		len.drawn.i <- length(relatedOntoID.direct.i)
		p.value.i <- NULL
		evdInfo.i <- NULL
		useCts.i <- NULL
		offsprings.i <- Onto2offspring[relatedOntoID.i]
		for (k in seq(relatedOntoID.i)) {
			relatedOntoID.ik <- relatedOntoID.i[k]
			
			## Determine whether to include the direct evidence of the terms themselves
			if (includeTestOnto) {
				offspringEvd.ik <- relatedOntoID.direct.i[relatedOntoID.direct.i %in% c(offsprings.i[[k]], relatedOntoID.ik)]
			} else {
				offspringEvd.ik <- relatedOntoID.direct.i[relatedOntoID.direct.i %in% offsprings.i[[k]]]
			}
			len.Evd.ik <- length(offspringEvd.ik)
			num.Onto2allEvd.ik <- num.Onto2allEvd[relatedOntoID.ik]

			useCts.i <- c(useCts.i, len.Evd.ik)
			if (len.Evd.ik >= minNumTh) {
				# p.value.ik <- phyper(len.Evd.ik, num.Onto2allEvd.ik, numAllEvd - num.Onto2allEvd.ik, len.drawn.i, lower.tail = FALSE)
				p.value.ik <- phyper(len.Evd.ik - 1, num.Onto2allEvd.ik, numAllEvd - num.Onto2allEvd.ik, len.drawn.i, lower.tail = FALSE)
			} else {
				p.value.ik <- 1
			}
			p.value.i <- c(p.value.i, p.value.ik)
			evdInfo.i <- c(evdInfo.i, list(offspringEvd.ik))
		}
		names(p.value.i) <- names(evdInfo.i) <- relatedOntoID.i

		## before sorting, we need to break the tie based on whether there are direct mapping or not
		tmp <- sort(table(relatedOntoID.direct.i), decreasing=TRUE)
		p.value.i <- c(p.value.i[names(tmp)], p.value.i[!(relatedOntoID.i %in% relatedOntoID.direct.i)])
		p.value.i <- sort(p.value.i, decreasing=FALSE)
		evdInfo.i <- evdInfo.i[names(p.value.i)]
		p.adjust.i <- p.adjust(p.value.i, method=fdr.adjust)
		p.value.sig.i <- p.value.i[(p.adjust.i < p.value.th)]
		ontoID.sig.i <- names(p.value.sig.i)

		if (length(ontoID.sig.i) > 0) {
			evdInfo.i <- evdInfo.i[names(p.value.sig.i)]
			p.adjust.i <- p.adjust.i[names(p.value.sig.i)]
			bestOntoID.i <- ontoID.sig.i[which.min(p.value.sig.i)]
			bestPvalue.i <- p.value.i[bestOntoID.i]
			bestP.adjust.i <- p.adjust.i[bestOntoID.i]
			
			names(p.value.sig.i) <- NULL
			names(p.adjust.i) <- NULL
			if (length(p.value.sig.i) > 1) {
				testSigInfo.i <- lapply(seq(ontoID.sig.i), function(ind) list(pValue=p.value.sig.i[ind], p.adjust=p.adjust.i[ind], evidence=evdInfo.i[[ind]]))
			} else {
				testSigInfo.i <- list(list(pValue=p.value.sig.i, p.adjust=p.value.sig.i, evidence=evdInfo.i[[1]]))
			}
			allSigOntoID.i <- c(allSigOntoID.i, ontoID.sig.i)
		}
		
		names(testSigInfo.i) <- allSigOntoID.i
		if (length(bestPvalue.i) > 0) {
			names(bestPvalue.i) <- names(bestP.adjust.i) <- bestOntoID.i
			bestOntoInfo.i <- list(pValue=bestPvalue.i, p.adjust=bestP.adjust.i)
		} else {
			bestOntoInfo.i <- NULL
		}
		
		testSigInfo.i <- list(list(allEvidence=gene2Onto[[gene.i]], sigOntoInfo=testSigInfo.i, bestOntoInfo=bestOntoInfo.i))
		testSigInfo <- c(testSigInfo, testSigInfo.i)		
	}
	names(testSigInfo) <- genes
	class(testSigInfo) <- 'geneFunSummarizeTest'
	attr(testSigInfo, 'parameters') <- list(rmOntoID=rmOntoID, p.value.th=p.value.th, fdr.adjust=fdr.adjust, minNumTh=minNumTh, includeTestOnto=includeTestOnto, directMapConstraint=directMapConstraint)
	return(testSigInfo)
}


## Simplify the significant ontology terms to a mini-set, which includes the non-overlapping most significant terms
simplifyGeneFunSummary <- function(geneFunSummarizeTestInfo, Onto.graph.closure, allOntoID.direct=NULL, fdr.adjust='none', p.value.th=10^-5) {
	
	genes <- names(geneFunSummarizeTestInfo)
	simplifyInfo <- NULL
	for (i in seq(genes)) {
		gene.i <- genes[i]
		geneFunSummarizeTestInfo.i <- geneFunSummarizeTestInfo[[i]]
		if (is.null(geneFunSummarizeTestInfo.i$sigOntoInfo)) {
			keptEvidences.i <- geneFunSummarizeTestInfo.i$allEvidence
			scores.i <- rep(1, length(keptEvidences.i))
			names(scores.i) <- keptEvidences.i
			simplifyInfo.i <- list(keptSigOntoID=NULL, keptEvidences=keptEvidences.i, scores=scores.i)
			simplifyInfo <- c(simplifyInfo, list(simplifyInfo.i))
			next
		}
		if (fdr.adjust == 'none') {
			pValues.i <- sapply(geneFunSummarizeTestInfo.i$sigOntoInfo, function(x) x$pValue)
		} else {
			pValues.i <- sapply(geneFunSummarizeTestInfo.i$sigOntoInfo, function(x) x$p.adjust)
		}
		sigPValues.i <- pValues.i[pValues.i < p.value.th]
		relatedOntoID.i <- names(sigPValues.i)
		if (!is.null(allOntoID.direct)) relatedOntoID.i <- relatedOntoID.i[relatedOntoID.i %in% allOntoID.direct]
		relatedEvidences.i <- lapply(geneFunSummarizeTestInfo.i$sigOntoInfo[relatedOntoID.i], function(x) x$evidence)
		allEvidences.i <- unique(geneFunSummarizeTestInfo.i$allEvidence)
		keptEvidences.i <- allEvidences.i[!(allEvidences.i %in% unlist(relatedEvidences.i))]
		keptSigOntoID.i <- NULL
		usedEvidences.i <- NULL
		while (length(relatedOntoID.i) > 0) {
			keptSigOntoID.ii <- relatedOntoID.i[which.min(sigPValues.i)]
			keptSigOntoID.i <- c(keptSigOntoID.i, keptSigOntoID.ii)
			relatedOntoID.i <- relatedOntoID.i[-which.min(sigPValues.i)]
			offspring.sig.i <- adj(Onto.graph.closure, keptSigOntoID.ii)[[1]]
			evidence.sig.ii <- geneFunSummarizeTestInfo.i$sigOntoInfo[[keptSigOntoID.ii]]$evidence
			usedEvidences.i <- unique(c(usedEvidences.i, evidence.sig.ii))
			
			## remove all offspring
			rmOntoInd.i <- which(relatedOntoID.i %in% offspring.sig.i)
			rmOntoID.i <- NULL
			if (length(rmOntoInd.i) > 0) {
				rmOntoID.i <- relatedOntoID.i[rmOntoInd.i]
				relatedOntoID.i <- relatedOntoID.i[-rmOntoInd.i]
			}
			if (length(relatedOntoID.i) == 0) break
			
			## remove all ancestors if they have the same evidences (direct mappings)
			for (j in seq(relatedOntoID.i)) {
				relatedOntoID.ij <- relatedOntoID.i[j]
				directOntoID.sig.ij <- relatedEvidences.i[[relatedOntoID.ij]]
				
				# if (all(directOntoID.sig.ij %in% evidence.sig.ii)) {
				if (all(directOntoID.sig.ij %in% usedEvidences.i)) {
					rmOntoID.i <- c(rmOntoID.i, relatedOntoID.ij)
				}
			}
			if (length(rmOntoID.i) > 0) {
				relatedOntoID.i <- relatedOntoID.i[!(relatedOntoID.i %in% rmOntoID.i)]
			}
		}
		## check each keptEvidences, whether they are ancestors of sigOntoIDs
		rmInd.i <- NULL
		for (j in seq(keptEvidences.i)) {
			offspring.ij <- adj(Onto.graph.closure, keptEvidences.i[j])[[1]]
			if (any (keptSigOntoID.i %in% offspring.ij)) {
				rmInd.i <- c(rmInd.i, j)
			}
		}
		if (!is.null(rmInd.i)) keptEvidences.i <- keptEvidences.i[-rmInd.i]
		scores.i <- round(-log10(pValues.i[c(keptSigOntoID.i, keptEvidences.i)]), 1)
		scores.i[scores.i < 1 | is.na(scores.i)] <- 1
		names(scores.i) <- c(keptSigOntoID.i, keptEvidences.i)
		simplifyInfo.i <- list(keptSigOntoID=keptSigOntoID.i, keptEvidences=keptEvidences.i, scores=scores.i)
		simplifyInfo <- c(simplifyInfo, list(simplifyInfo.i))
	}
	names(simplifyInfo) <- genes
	attr(simplifyInfo, 'pValueT') <- p.value.th
	return(simplifyInfo)
}


saveGeneFunSummary <- function(geneFunSummarizeResult, simplifyInfo=NULL, addFDR=FALSE, species='human', ID2Name=NULL, fileName='geneSummarization.xls') {

	genes <- names(geneFunSummarizeResult)
	allEvidence <- sapply(geneFunSummarizeResult, function(x) paste(x$allEvidence, collapse='; '))
	bestOntoInfo <- sapply(geneFunSummarizeResult, function(x) {
		if (length(x$bestOntoInfo) == 0) {
			return('')
		} else {
			bestOntoInfo.x <- x$bestOntoInfo
			bestOntoID.x <- names(bestOntoInfo.x[[1]])
			pValue.x <- bestOntoInfo.x$pValue
			pAdjust.x <- bestOntoInfo.x$p.adjust
			if (is.null(ID2Name)) {
				if (addFDR) {
					return(paste(paste(bestOntoID.x, ' (', signif(-log10(pValue.x),3), ', ', signif(-log10(pAdjust.x),3), ')', sep=''), collapse='; '))
				} else {
					return(paste(paste(bestOntoID.x, ' (', signif(-log10(pValue.x),3), ')', sep=''), collapse='; '))
				}
			} else {
				bestOntoName.x <- ID2Name[bestOntoID.x]
				if (addFDR) {
					return(paste(paste(bestOntoName.x, ' (', bestOntoID.x, ')', ' (', signif(-log10(pValue.x),3), ', ', signif(pAdjust.x,3),')', sep=''), collapse='; '))
				} else {
					return(paste(paste(bestOntoName.x, ' (', bestOntoID.x, ')', ' (', signif(-log10(pValue.x),3), ')', sep=''), collapse='; '))
				}
			}
		}
	})
	if (!is.null(simplifyInfo)) {
		miniSetInfo <- sapply(seq(simplifyInfo), function(i) {
			ans.i <- geneFunSummarizeResult[[i]]
			sim.i <- simplifyInfo[[i]]
			simOnto.i <- sim.i$keptSigOntoID
			if (!is.null(simOnto.i)) {
				simPvalue.i <- sapply(ans.i$sigOntoInfo[simOnto.i], function(x) x$pValue)
				if (addFDR) {
					simFDR.i <- sapply(ans.i$sigOntoInfo[simOnto.i], function(x) x$p.adjust)
					if (is.null(ID2Name)) {
						sigInfo.i <- paste(paste(simOnto.i, " (", signif(-log10(simPvalue.i),3), ", ", signif(-log10(simFDR.i),3), ")", sep=""), collapse="; ")
					} else {
						sigInfo.i <- paste(paste(ID2Name[simOnto.i], " (", signif(-log10(simPvalue.i),3), ", ", signif(-log10(simFDR.i),3), ")", sep=""), collapse="; ")
					}
				} else {
					sigInfo.i <- paste(paste(ID2Name[simOnto.i], " (", signif(-log10(simPvalue.i) ,3), ")", sep=""), collapse="; ")
				}
			} else {
				sigInfo.i <- ""
			}
			keptOnto.i <- sim.i$keptEvidences
			if (!is.null(keptOnto.i)) {
				if (is.null(ID2Name)) {
					keptOntoInfo.i <- paste(paste(keptOnto.i, " (1)", sep=""), collapse="; ")
				} else {
					keptOntoInfo.i <- paste(paste(ID2Name[keptOnto.i], " (1)", sep=""), collapse="; ")
				}
				if (sigInfo.i != "") {
					sigInfo.i <- paste(sigInfo.i, keptOntoInfo.i, sep="; ")					
				} else {
					sigInfo.i <- keptOntoInfo.i				
				}				
			} 
			return(sigInfo.i)
		})		
	}
	summaryInfo <- sapply(geneFunSummarizeResult, function(x) {
		if (length(x$sigOntoInfo) == 0) {
			return('')
		} else {
			sigOntoID.x <- names(x$sigOntoInfo)
			pValue.x <- sapply(x$sigOntoInfo, function(x) x$pValue)
			pAdjust.x <- sapply(x$sigOntoInfo, function(x) x$p.adjust)
			return(paste(paste(sigOntoID.x, ' (', signif(-log10(pValue.x),3), ', ', signif(-log10(pAdjust.x),3),')', sep=''), collapse='; '))
		}
	})
	lib <- switch(species,
		'rat'='org.Rn.eg.db',
		'human'='org.Hs.eg.db',
		'mouse'='org.Mm.eg.db',
		'yeast'='org.Sc.eg.db',
		'fly'='org.Dm.eg.db'
		 )
	if (require(lib, character=TRUE)) {
		geneSymbol <- getSYMBOL(genes, lib)
		if (is.null(simplifyInfo)) {
			geneSummary <- data.frame(geneID=genes, geneSymbol=geneSymbol, bestOntology=bestOntoInfo, enrichedOntology=summaryInfo, allEvidence=allEvidence)
		} else {
			geneSummary <- data.frame(geneID=genes, geneSymbol=geneSymbol, bestOntology=bestOntoInfo, miniSetOntology=miniSetInfo, enrichedOntology=summaryInfo, allEvidence=allEvidence)
		}
	} else {
		if (is.null(simplifyInfo)) {
			geneSummary <- data.frame(geneID=genes, bestOntology=bestOntoInfo, enrichedOntology=summaryInfo, allEvidence=allEvidence)
		} else {
			geneSummary <- data.frame(geneID=genes, bestOntology=bestOntoInfo, miniSetOntology=miniSetInfo, enrichedOntology=summaryInfo, allEvidence=allEvidence)
		}
	}
	if (!is.null(fileName)) {
		write.table(geneSummary, file=fileName, sep='\t', row.names=FALSE, quote=FALSE)
	}
	return(invisible(geneSummary))
}


plotGeneFunSummary <- function(geneFunSummarizeResult, onto.graph, selGene=NULL, onto.graph.closure=NULL, allOntoID.direct=NULL, fdr.adjust='none', showMiniSet=TRUE,  highlightBest=TRUE,  miniSetPvalue=10^-5, ID2Name=NULL, savePrefix='', geneSymbol=TRUE, saveImage=TRUE, selectionMethod=c('simple', 'bestOnly',  'all'), lib='org.Hs.eg.db', ...) {

	selectionMethod <- match.arg(selectionMethod)
	genes <- names(geneFunSummarizeResult)
	if (!is.null(selGene)) genes <- intersect(genes, selGene)
	geneGraph <- NULL
	sigGene <- NULL
	for (i in seq(genes)) {
		gene.i <- genes[i]
		geneFunSummarize.i <- geneFunSummarizeResult[[gene.i]]
		if (fdr.adjust == 'none') {
			pValue.i <- sapply(geneFunSummarize.i$sigOntoInfo, function(x) x$pValue)
		} else {
			pValue.i <- sapply(geneFunSummarize.i$sigOntoInfo, function(x) x$p.adjust)
		}
		# pValue.i <- pValue.i[pValue.i <= p.value.th]
		if (length(pValue.i) == 0) next
		sigGene <- c(sigGene, gene.i)
		relatedOntoID.i <- geneFunSummarize.i$allEvidence
		if (showMiniSet) {
			simplifyInfo <- simplifyGeneFunSummary(geneFunSummarizeResult[i], onto.graph.closure, allOntoID.direct=allOntoID.direct, p.value.th=miniSetPvalue)
			bestOntoID.i <- names(simplifyInfo[[1]]$scores)
		} else if (highlightBest) {
			bestOntoID.i <- names(geneFunSummarize.i$bestOntoInfo[[1]])			
		} else {
			bestOntoID.i <- NULL
		}
		if (geneSymbol && require(lib, character=TRUE)) gene.i <- getSYMBOL(gene.i, lib)
		if (saveImage) {
			if (is.null(savePrefix) || savePrefix == '') {			
				saveImageName.i <- paste('ontologyGraph', gene.i, sep="_")
			} else {
				saveImageName.i <- paste(savePrefix, gene.i, sep="_")
			}
		} else {
			saveImageName.i <- NULL
		}

		geneGraph.i <- plotOntologyGraph(pValue.i, relatedOntoID=relatedOntoID.i, onto.graph=onto.graph, bestOntoID=bestOntoID.i, onto.graph.closure=onto.graph.closure, ID2Name=ID2Name, saveImageName=saveImageName.i, selectionMethod=selectionMethod, ...)
		geneGraph <- c(geneGraph, list(geneGraph.i))
	}
	names(geneGraph) <- sigGene
	return(invisible(geneGraph))
}


plotOntologyGraph <- function(onto.pValue, relatedOntoID, onto.graph, bestOntoID=NULL, onto.graph.closure=NULL, rootID="DOID:4", ID2Name=NULL, p.value.th=0.01, fillColor='white', colorLevel=seq(2,20,by=2), relative.color=TRUE, fontsize=15, colorMap=colorRampPalette(c('white', 'red'))(length(colorLevel)+1), selectionMethod=c('simple', 'bestOnly', 'all'), omitNode=TRUE, saveImageName=NULL) {
	
	selectionMethod <- match.arg(selectionMethod)
	if (selectionMethod == 'all') omitNode <- FALSE
	if (is.null(onto.graph.closure)) onto.graph.closure <- transitive.closure(onto.graph)

	r.onto.graph.closure <- reverseEdgeDirections(onto.graph.closure)
	p.value <- onto.pValue[onto.pValue <= p.value.th]
	
	if (length(p.value) == 0) {
		cat('No significant nodes were selected!\n')
		sigID <- NULL
	} else {
		sigID <- names(p.value)
		if (is.null(sigID)) stop('onto.pValue should include the names (IDs) of the node!\n')
	}
	
	if (!is.null(rootID) && selectionMethod == 'all') {
		relatedID.sig <- adj(onto.graph.closure, c(rootID, sigID))
	} else {
		relatedID.sig <- adj(onto.graph.closure, sigID)	
	}
	relatedID.sig <- unique(unlist(relatedID.sig))
	relatedID.r <- adj(r.onto.graph.closure, relatedOntoID)
	relatedID.r <- unique(unlist(relatedID.r))
	commID <- unique(c(intersect(relatedID.sig, relatedID.r), relatedOntoID, sigID, bestOntoID))
	if (length(commID) < 3)  return(NULL)

	subg <- subGraph(commID, onto.graph)
	
	if (selectionMethod %in% c('simple', 'bestOnly')) {
		subg.r <- reverseEdgeDirections(subg)
		
		if (selectionMethod == 'simple') {
			if (!is.null(bestOntoID)) {
				acc.best <- acc(subg, bestOntoID)
				sigID <- sigID[sigID %in% c(bestOntoID, unlist(lapply(acc.best, names)))]
				interestedID <- c(sigID, bestOntoID)
			} else {
				interestedID <- sigID
			}
		} else {
			interestedID <- bestOntoID
		}
		keptID <- unique(c(relatedOntoID, interestedID))
		selInd <- interestedID %in% nodes(subg)
		interestedID <- interestedID[selInd]
		## add code here:

		addEdges <- NULL
		if (omitNode) {
			subg.final <- subGraph(keptID, subg)
			acc.all <- acc(subg, keptID)			
			for(i in seq(acc.all)) {
				acc.i <- acc.all[[i]]
				from.i <- keptID[i]
				acc.IDs <- names(acc.i)
				# remove the nodes whihc are the offspring of keptID				
				acc.IDs <- acc.IDs[acc.IDs %in% keptID]
				# selID.i <- acc.IDs[acc.IDs %in% sigID]
				if (length(acc.IDs) > 0) {	
					acc.IDs <- acc.IDs[!(acc.IDs %in% unlist(lapply(acc.all[acc.IDs], names)))]					
				} 
				acc.IDs <- acc.IDs[acc.i[acc.IDs] > 1]
				if (length(acc.IDs) > 0) {
					# add edges with indirect links
					for (acc.ID.i in acc.IDs) {
						subg.final <- addEdge(from.i, acc.ID.i, subg.final, acc.i[acc.ID.i]) # weights is the steps between nodes
						addEdges <- c(addEdges, paste(from.i, acc.ID.i, sep="~"))
					}
				}
			}
			
		} else {
			interestedID <- unique(c(interestedID, relatedOntoID))
			
			acc.all <- acc(subg, interestedID)
			closest.dist <- lapply(acc.all, function(x) {
				acc.IDs <- names(x)
				selID.i <- acc.IDs[acc.IDs %in% interestedID]
				if (length(selID.i) == 0) {
					temp = 0
					names(temp) = 'NA'
					return(temp)
				} else {
					return(x[selID.i][which.min(x[selID.i])])
				}
			})
			closest.ID <- sapply(closest.dist, names)
			closest.dist <- unlist(closest.dist)
			keptID <- c(keptID, closest.ID[closest.dist == 1])
			closestID <- closest.ID[closest.dist > 1]
			if (length(closestID) > 0) {
				sp <- sp.between(subg, names(closestID), closestID)
				pathNodes <- unlist(sapply(sp, function(x) x$path_detail))
				keptID <- unique(c(keptID, pathNodes))
			}
			# deal with the reverse direction
			acc.all.r <- acc(subg.r, interestedID)
			closest.dist.r <- lapply(acc.all.r, function(x) {
				acc.IDs <- names(x)
				selID.i <- acc.IDs[acc.IDs %in% interestedID]
				if (length(selID.i) == 0) {
					temp = 0
					names(temp) = 'NA'
					return(temp)
				} else {
					return(x[selID.i][which.min(x[selID.i])])
				}
			})
			closest.ID.r <- sapply(closest.dist.r, names)
			closest.dist.r <- unlist(closest.dist.r)
			keptID <- c(keptID, closest.ID.r[closest.dist.r == 1])
			closestID.r <- closest.ID.r[closest.dist.r > 1]
			if (length(closestID.r) > 0) {
				sp <- sp.between(subg.r, names(closestID.r), closestID.r)
				pathNodes.r <- unlist(sapply(sp, function(x) x$path_detail))
				keptID <- c(keptID, pathNodes.r)
			}
			if (!is.null(bestOntoID)) {
				relatedID.best <- adj(onto.graph.closure, bestOntoID)	
				relatedID.best <- unique(unlist(relatedID.best))
				commID.best <- c(intersect(relatedID.best, relatedID.r), bestOntoID)
				keptID <- c(keptID, commID.best)
			}
			subg.final <- subGraph(unique(keptID), onto.graph)			
		}	

	} else {
		subg.final <- subg
	}

	## -----------------------
	## plot the figure with the significant ID highlighted
	# colorMap <- colorRampPalette(c('white', 'red'))(32)
	nn <- nodes(subg.final)
	shape <- rep('ellipse', length(nn))
	names(shape) <- nn
	if (length(onto.pValue) > 0) {
		log.p <- -log10(onto.pValue)
		best.p.value <- round(log.p)
		if (relative.color) {
			# log.p <- (log.p - -log10(p.value.th))/(max(log.p) - -log10(p.value.th))
			log.p <- (log.p - 2)/(max(log.p) - 2)
			color.ind <- round(log.p * (length(colorLevel) - 1)) + 2
		} else {
			color.ind <- as.numeric(cut(log.p, c(colorLevel, Inf)))	+ 1		
		}
		selNode.color <- colorMap[color.ind]
		names(selNode.color) <- names(onto.pValue)
		if (!is.null(bestOntoID)) {
			#bestOntoID <- names(onto.pValue)[which.min(onto.pValue)]
			shape[bestOntoID] <- 'rectangle'
		}
	} else {
		best.p.value <- 0
	}
	fillColor <- rep(fillColor, length(relatedOntoID))
	names(fillColor) <- relatedOntoID
	fillColor <- c(selNode.color, fillColor)

	directNodeCol <- rep('green', length(relatedOntoID))
	directNodeLW <- rep(1.5, length(relatedOntoID))
	directNodeLty <- rep(1, length(relatedOntoID))
	names(directNodeCol) <- names(directNodeLW) <- names(directNodeLty) <- relatedOntoID
	lwd <- rep(1, length(nodes(subg.final)))
	names(lwd) <- nodes(subg.final)
	lwd[relatedOntoID] <- 1.5

	layoutAttrs <- list(node=list(shape=shape, fixedsize=FALSE))
	if (!is.null(ID2Name)) {
		layoutNodeAttrs <- list(label=ID2Name[nodes(subg.final)])
	} else {
		layoutNodeAttrs <- list(label=nodes(subg.final))
	}
	g1 <- layoutGraph(subg.final, attrs=layoutAttrs, nodeAttrs=layoutNodeAttrs)
	if (length(unlist(edges(g1))) == 0) return(invisible(list(graph=subg.final,layout=g1)))
	if (!is.null(saveImageName)) {
		if (relative.color) saveImageName <- paste('relativeColor.', saveImageName, sep='')
		saveImageName <- sub('\\.pdf$', '', saveImageName)
		saveImageName <- paste(saveImageName, "_", best.p.value, '.pdf', sep='')
		pdf(file=saveImageName, width=8, height=4)		
	}
	nodeRenderAttrs <- list(fill=fillColor, lwd=lwd, col=directNodeCol, lty=directNodeLty, fontsize=fontsize, shape=shape)
	nodeRenderInfo(g1) <- nodeRenderAttrs

	edgeRenderAttrs <- NULL
	if (!is.null(addEdges)) {
		lty <- rep('dashed', length(addEdges))
		names(lty) <- addEdges
		edgeRenderAttrs <- list(lty=lty)
		edgeRenderInfo(g1) <- edgeRenderAttrs
	}
	renderGraph(g1)
	if (!is.null(saveImageName)) dev.off()

	return(invisible(list(graph=subg.final, layoutAttrs=layoutAttrs, layoutNodeAttrs=layoutNodeAttrs, layoutEdgeAttrs=NULL, nodeRenderAttrs=nodeRenderAttrs, edgeRenderAttrs=edgeRenderAttrs)))
}


plotGraph <- function(graph, layoutAttrs=NULL, layoutNodeAttrs=NULL, layoutEdgeAttrs=NULL, nodeRenderAttrs=NULL, edgeRenderAttrs=NULL) {

	g1 <- layoutGraph(graph, attrs=layoutAttrs, nodeAttrs=layoutNodeAttrs, edgeAttrs=layoutEdgeAttrs)
	if (!is.null(nodeRenderAttrs))  nodeRenderInfo(g1) <- nodeRenderAttrs
	if (!is.null(edgeRenderAttrs))  edgeRenderInfo(g1) <- edgeRenderAttrs
	renderGraph(g1)

	return(invisible(list(graph=graph,layoutAttrs=layoutAttrs, layoutNodeAttrs=layoutNodeAttrs, layoutEdgeAttrs=layoutNodeAttrs, nodeRenderAttrs=nodeRenderAttrs, edgeRenderAttrs=edgeRenderAttrs)))	
}
