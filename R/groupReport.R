`groupReport` <- 
function(dataMatrix, gAList, topCat=10, methodOfCluster=c('mds', 'sort'), matrixOfHeatmap=NULL, clusterTable=c('geneNum', 'pvalue', NULL), catTerm=TRUE, 
		fileName = "multiConceptsGenes.html", title='Multigroup Genes Concepts Analysis', catType=c('GO', 'KEGG', 'DOLITE', 'REACTOME.PATH', 'CABIO.PATH', 'Unknown'), reverseOfCluster=FALSE,  colorValueColumn = NULL, 
		annLib=c('org.Hs.eg.db', 'org.Rn.eg.db', 'org.Mm.eg.db', 'org.Dm.eg.db'), nameLength=94, addID=TRUE, interactive=FALSE, bgColor='#ffffcc', keepCytoscapeFiles=TRUE, ...) { 
	catType <-match.arg(catType)
	annLib <- match.arg(annLib)
	
	showCats <- c(1:topCat)
	
	k <- 0
	while (TRUE) {
		if ((all(!file.exists(c(paste(fileName, '_', k, '.html', sep=''), paste(paste(fileName, '_', k, '.html', sep=''), 'files', sep='.'), paste(paste(fileName, '_', k, '.html', sep=''), 'cytoscapeWebFiles', sep='.')))) & (k > 0)) | 
			(all(!file.exists(c(fileName, paste(fileName, 'files', sep='.'), paste(fileName, 'cytoscapeWebFiles', sep='.')))) & (k == 0))) {
			if (k > 0) {
				fileName <- paste(fileName, '_', k, '.html', sep='')
				print(paste('New file is', fileName))
			}
			break;
		} else {
			if (k > 0) print(paste(paste(c(paste(fileName, '_', k, '.html', sep=''), paste(paste(fileName, '_', k, '.html', sep=''), 'files', sep='.'), paste(paste(fileName, '_', k, '.html', sep=''), 'cytoscapeWebFiles', sep='.')), 
						collapse=' and/or '), 'are in current location, they will not be overwritten.'))
			else print(paste(paste(c(fileName, paste(fileName, 'files', sep='.'), paste(fileName, 'cytoscapeWebFiles', sep='.')), collapse=' and/or '), 'are in the current location, they will not be overwritten.'))
		}
		k <- k + 1
	}
	
	if (!is.null(colorValueColumn) & (length(colorValueColumn) > 1) & (length(colorValueColumn) != length(gAList))) stop('Specified colorValueColumns can not match the given gAList! Aborting ...') 
	
	if (length(colorValueColumn) == 1) colorValueColumn <- rep(colorValueColumn, length(gAList))
	
	if (!is.null(colorValueColumn)) {
		allColNames <- lapply(gAList, function(x) return(colnames(getGeneInput(x))))
		validColNames <- c()
		for (i in 1:length(colorValueColumn)) validColNames <- c(validColNames, colorValueColumn[i] %in% allColNames[[i]])
		if (!(all(validColNames))) stop(paste("The given colorValueColumn(s)[position ", paste(which(!validColNames), sep="", collapse=" "), "] might not be the column names of geneInput slots of the given gAList! Aborting ...", sep=''))
	}
	
	if (is.matrix(dataMatrix) | is.data.frame(dataMatrix) | is.vector(dataMatrix)) {
		dataMatrix <- as.matrix(dataMatrix)
		if ((NA %in% rownames(dataMatrix)) | is.null(rownames(dataMatrix))) print('Warning: NA or NULL might be in rownames of dataMatrix!')
		if ((NA %in% colnames(dataMatrix)) | is.null(colnames(dataMatrix))) print('Warning: NA or NULL might be in colnames of dataMatrix!')
	} else stop('Input is not a valid matrix!')
	
	HTwrap <- function(x, tag = "TD", scripts=NULL) {
        if (is.null(scripts)) return(paste("<", tag, ">", x, "</", tag, ">", sep = ""))
		else return(paste("<", tag, ' ', scripts, ">", x, "</", tag, ">", sep = ""))
    }

    if (!is.null(clusterTable) & !is.null(matrixOfHeatmap)) {
		if (clusterTable == 'geneNum') {
			if (is.numeric(dataMatrix)) {
				baseM <- dataMatrix[1:(dim(dataMatrix)[1]-1),]
			}
			else stop('The input Matrix is not numeric!')
		} else {
			baseM <- matrixOfHeatmap
		}

		if (methodOfCluster == 'mds') {
			distBaseM <- as.matrix(dist(baseM, method = "euclidean", upper=TRUE))
			distBaseM[distBaseM < 1e-6] <- 1e-6
			diag(distBaseM) <- 0
			clusterResult <- isoMDS(distBaseM, k=1)$points
			clusterOrder <- rownames(clusterResult)[order(clusterResult, decreasing=reverseOfCluster)]
			originalIndexMatrix <- originalIndexMatrix[clusterOrder,]
			matrixOfHeatmap <- matrixOfHeatmap[clusterOrder,]
			tempM <- dataMatrix[c(clusterOrder, rownames(dataMatrix)[dim(dataMatrix)[1]]),]
			dataMatrix <- tempM
		} else {
			clusterOrder <- sort(rownames(baseM), ...)
			tempM <- dataMatrix[c(clusterOrder, rownames(dataMatrix)[dim(dataMatrix)[1]]),]
			dataMatrix <- tempM
		}
	}
 
    indexM <- which((is.numeric(dataMatrix) & (dataMatrix > 0)) | (is.character(dataMatrix) & (dataMatrix != "0 (1)")), arr.ind = TRUE)
    if (catType != 'Unknown') {
	 	indexNames <- rownames(indexM)[1:dim(indexM)[1]]
		indexTerms <- getCategoryTerms(indexNames[which(!indexNames %in% "Genes / Group")], catType = catType, missing = "keep")
		trueTerms <- paste('Genes in', names(indexTerms), '::', indexTerms, 'and Group', colnames(dataMatrix)[indexM[1:(dim(indexM)[1]),2][which(!(names(indexM[1:(dim(indexM)[1]),2]) %in% "Genes / Group"))]])
		names(trueTerms) <- names(indexTerms)
		pseudoTerms <- paste('Genes in Group', colnames(dataMatrix))
		names(pseudoTerms) <- rep('Genes / Group', length(pseudoTerms))
		finalTerms <- rep('', length(indexNames))
		finalTerms[which(!indexNames %in% "Genes / Group")] <- trueTerms
		finalTerms[which(indexNames %in% "Genes / Group")] <- pseudoTerms
		names(finalTerms) <- NULL
		tableNames <- c('Concepts-Genes Table', 'Concepts Relation', paste('Group', colnames(dataMatrix), 'Concepts-Gene network'), finalTerms)
		#gsub('Genes / Group :: NA and ', replacement='', paste('Genes in', rownames(indexM)[1:dim(indexM)[1]], '::', getCategoryTerms(rownames(indexM)[1:dim(indexM)[1]], catType=catType, missing='keep'), 'and Group', colnames(dataMatrix)[indexM[1:dim(indexM)[1],2]]))
	} else {
		tableNames <- c('Concepts-Genes Table', 'Concepts Relation', paste('Group', colnames(dataMatrix), 'Concepts-Gene network'), gsub('Genes / Group and ', replacement='', 
							paste('Genes in', rownames(indexM)[1:dim(indexM)[1]], 'and Group', colnames(dataMatrix)[indexM[1:dim(indexM)[1],2]])))
	}
    attr <- c(catType, rep('Pictures', (1+dim(dataMatrix)[2])), rep('Entrez', dim(indexM)[1]))
	outFile <- file(fileName, "w")
	cat('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN""http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">', file = outFile, sep = "\n")
    
	if (dir.create(paste(fileName, '.files', sep=''), showWarnings=FALSE) | file.exists(paste(fileName, '.files', sep=''))) {
		if (interactive) {
			if (file.copy(paste(system.file('External', package='GeneAnswers'), '/js', sep=''), paste(fileName, '.files', sep=''), recursive=TRUE)) print(paste('Copy js to', paste(fileName, '.files', sep='')))
			else { close(outFile); stop('js subdirectory can not be copied to ', paste(fileName, '.files', sep=''), '. Aborting ...') }
			if (file.copy(paste(system.file('External', package='GeneAnswers'), '/swf', sep=''), paste(fileName, '.files', sep=''), recursive=TRUE)) print(paste('Copy swf to', paste(fileName, '.files', sep='')))
			else { close(outFile); stop('swf subdirectory can not be copied to ', paste(fileName, '.files', sep=''), '. Aborting ...') }
			if (file.copy(paste(system.file('External', package='GeneAnswers'), '/NetworkTransformer.jar', sep=''), '.', overwrite=TRUE)) {
				Sys.chmod('NetworkTransformer.jar', mode='0777')
				print(paste('Copy NetworkTransformer.jar to', paste(fileName, '.files', sep='')))
			}
			else { close(outFile); stop('NetworkTransformer.jar subdirectory can not be copied to ', paste(fileName, '.files', sep=''), '. Aborting ...') }

			headInfo <- paste(HTwrap(title, tag = "title"), '\n\t', HTwrap('', tag='script', scripts=paste('type="text/javascript" src="', paste(fileName, '.files', sep=''), '/js/cytoscape_web/json2.min.js"', sep='')), '\n\t', 
								HTwrap('', tag='script', scripts=paste('type="text/javascript" src="', paste(fileName, '.files', sep=''), '/js/cytoscape_web/AC_OETags.min.js"', sep='')), '\n\t', 
								HTwrap('', tag='script', scripts=paste('type="text/javascript" src="', paste(fileName, '.files', sep=''), '/js/cytoscape_web/cytoscapeweb.min.js"', sep='')), '\n', sep='')
		    headInfo <- paste(headInfo, paste('\t', HTwrap('', tag='script', scripts=paste('type="text/javascript" src="', paste(fileName, '.files', sep=''), '/edges_', c(1:length(gAList)), '.js"', sep='')),'\n', sep='', collapse=''), 
							'\t\t<script type="text/javascript">\n\t\twindow.onload=function() {\n', sep='')
			headInfo <- paste(headInfo, paste('\t\t\tvar div_id', c(1:length(gAList)), ' = "cytoscapeweb', c(1:length(gAList)), '";\n', sep='', collapse=''), 
							'\t\t\tvar options = {\n\t\t\t\tswfPath: "', paste(fileName, '.files', sep=''), '/swf/CytoscapeWeb",\n\t\t\t\tflashInstallerPath: "', paste(fileName, '.files', sep=''), '/swf/playerProductInstall"\n\t\t\t};\n', sep='')
			headInfo <- paste(headInfo, paste('\t\t\tvar vis', c(1:length(gAList)), ' = new org.cytoscapeweb.Visualization(div_id', c(1:length(gAList)), ', options);\n\t\t\tvis', c(1:length(gAList)), '.draw({ network: xml',c(1:length(gAList)), ' });\n', sep='', collapse=''), 
							'\t\t};</script>\n\t', HTwrap('   html, body { height: 100%; width: 100%; padding: 0; margin: 10; }\n\t#cytoscapeweb { width: 100%; height: 100%; }', tag='style'), sep='')
		    cat("<html  xmlns='http://www.w3.org/1999/xhtml' xml:lang='en' lang='en'>", HTwrap(headInfo, tag = "head"), '<body bgcolor="', bgColor, '">', file = outFile, sep = "\n")
		} else {
			cat("<html  xmlns='http://www.w3.org/1999/xhtml' xml:lang='en' lang='en'>", HTwrap(HTwrap(title, tag = "title"), tag = "head"), '<body bgcolor="', bgColor, '">', file = outFile, sep = "\n")
		}
	} else {
		close(outFile)
		stop(paste('Failure to create', paste(fileName, '.files', sep=''), 'subdirectory! Aborting ...'))
	}
	

    if (!missing(title)) cat("<center><H1 align=\"center\">", title, " </H1></center>\n", file = outFile, sep = "\n")
	cat('</body>', "</html>", sep = "\n", file = outFile)
	
	cat('<H1><font face="courier" size="4"><a name="h-0">Multigroup Genes Concepts Analysis</a></font></H1>\n<font face="courier" size="1">Generated by Bioconductor package <a HREF="http://www.bioconductor.org/packages/release/bioc/html/GeneAnswers.html">GeneAnswers</a> at ', format(Sys.Date(), "%a, %b %d %Y"), '</font>\n<p>  </p>',file = outFile, sep = "")
	
	for (i in 1:(2+dim(dataMatrix)[2])) {
		cat(paste('<li><font face="courier" size="2"><a href="', paste('#h-', i, sep=''), '">', sep=''), tableNames[i], '</a></font>\n<p></p>\n',file = outFile, sep = "")
	}
	for (i in 1:length(tableNames)) {
		if (attr[i] != 'Pictures') {
			if (attr[i] == 'Entrez') {
				a <- i - 2 - dim(dataMatrix)[2] 
				tempGeneInput <- getGeneInput(gAList[[indexM[a,2]]])
				if (indexM[a,1] == dim(dataMatrix)[1]){
					# generate all genes table
					temp <- as.matrix(tempGeneInput)
					colnames(temp) <- colnames(tempGeneInput)
					rownames(temp) <- temp[,1]
					if (dim(temp)[2] > 1) {
	                    if (dim(temp)[1] > 1) {
							tempColNames <- colnames(temp)
							temp <- temp[!duplicated(temp),]
							if (!is.matrix(temp)) {
								temp <- matrix(temp, ncol=length(tempColNames))
								rownames(temp) <- temp[,1]
								colnames(temp) <- tempColNames
							}
						}
						if (dim(temp)[2] > 1) {
							tempColNames <- colnames(temp)[2:dim(temp)[2]]
							tempRowNames <- rownames(temp)
							temp <- as.matrix(temp[,2:dim(temp)[2]])
							colnames(temp) <- tempColNames
							rownames(temp) <- tempRowNames
						}
					}
  					.drawHTMLtable(temp, outFile, tableName=tableNames[i], tableLink=paste('h-', i, sep=''), catType=attr[i], species=getAnnLib(gAList[[1]]), lastRowLink=TRUE, highlightLastRow=FALSE, topCat=0, IDCols=2)
				} else {
					#generate genes in concept table
					tempGenes <- unlist(getGenesInCategory(gAList[[indexM[a,2]]])[rownames(dataMatrix)[indexM[a,1]]])
					
					temp <- as.matrix(tempGeneInput[tempGeneInput[,1] %in% tempGenes,])
					rownames(temp) <- temp[,1]
					colnames(temp) <- colnames(tempGeneInput)
					if (dim(temp)[2] > 1) {
						if (dim(temp)[1] > 1) {
							tempColNames <- colnames(temp)
							temp <- temp[!duplicated(temp),]
							if (!is.matrix(temp)) {
								temp <- matrix(temp, ncol=length(tempColNames))
								rownames(temp) <- temp[,1]
								colnames(temp) <- tempColNames
							}
						}
						if (dim(temp)[2] > 1) {
							tempColNames <- colnames(temp)[2:dim(temp)[2]]
							tempRowNames <- rownames(temp)
							temp <- matrix(temp[,2:dim(temp)[2]], ncol=length(tempColNames))
							colnames(temp) <- tempColNames
							rownames(temp) <- tempRowNames
						}
					}
					#### generate a hyperlink by table name based on catType
					.drawHTMLtable(temp, outFile, tableName=tableNames[i], tableLink=paste('h-', i, sep=''), catType=attr[i], species=getAnnLib(gAList[[1]]), lastRowLink=TRUE, highlightLastRow=FALSE, topCat=0, IDCols=2)
				}
            } else {
				.drawHTMLtable(dataMatrix, outFile, tableName=tableNames[i], tableLink=paste('h-', i, sep=''), catType=attr[i], species=annLib, matrixOfHeatmap=matrixOfHeatmap, topCat=max(showCats), displayText=TRUE)
			} 
		} else {
			if (i == 2) {
				if (catType == 'GO') tempConceptsRelationFileName <- paste(tableNames[i], '_', getCategoryType(gAList[[i-1]]), '.png', sep='')
				else tempConceptsRelationFileName <- paste(tableNames[i], '_', catType, '.png', sep='')
				setwd(paste(fileName, '.files', sep=''))
				png(filename=tempConceptsRelationFileName, width=1000, height=1000, bg='#ffffcc')
				if (is.matrix(dataMatrix[c(1:(dim(dataMatrix)[1]-1)),]))  {
					.catsCluster(dataMatrix[c(1:(dim(dataMatrix)[1]-1)),], gAList, catType=catType, nameLength=nameLength, addID=addID)
				} else {
					newTempM <- dataMatrix[c(1:(dim(dataMatrix)[1]-1)),]
					.catsCluster(matrix(newTempM, ncol=1, dimnames=list(c(names(newTempM)), c(names(gAList)))), gAList, catType=catType, nameLength=nameLength, addID=addID)
				}
				dev.off()
				setwd('..')
				cat(paste('<H2 align=center><font face="courier" size="2"><A name="', paste('h-', i, sep=''), '">', tableNames[i], "</a></font></H2>", sep=''), file = outFile, sep = "\n") 
				cat('<center><IMG src="', paste(fileName, '.files', sep=''), '/', tempConceptsRelationFileName, '"></center>\n', file = outFile, sep="")
			} else {
				if (!is.null(matrixOfHeatmap) & is.numeric(showCats)) {
					temp <- matrixOfHeatmap[,(i-2)]
					drawCats <- names(sort(temp))[showCats]
				} else drawCats <- showCats
				if (is.character(drawCats)) drawCats <- drawCats[drawCats %in% rownames(getEnrichmentInfo(gAList[[i-2]]))]
				else drawCats <- rownames(getEnrichmentInfo(gAList[[i-2]]))[intersect(showCats, c(1:dim(getEnrichmentInfo(gAList[[i-2]]))[1]))]
				if (!is.null(drawCats)) {
					##if (catType == 'GO') tempConceptFileName <- paste(tableNames[i], '_', getCategoryType(gAList[[i-2]]), '.png', sep='') 
					##else  tempConceptFileName <- paste(tableNames[i], '_', catType, '.png', sep='')
					tempCrossTableFileName <- paste(tableNames[i], '_', getCategoryType(gAList[[i-2]]), '_CrossTable.png', sep='')
					if (interactive) {
						g <- geneAnswersConceptNet(gAList[[i-2]], centroidSize='pvalue', colorValueColumn = colorValueColumn[i-2], output='none', showCats=drawCats, catTerm=catTerm, geneSymbol=TRUE, bgColor=bgColor, symmetry=FALSE)
						graphInfo <- c('vertex.attributes'=list(as.data.frame(cbind('NODES'=V(g)$label, 'NODE_FILL_COLOR'=V(g)$color, 'NODE_SIZE'=3*V(g)$size, 
																'NODE_LABEL_COLOR'=rep('#666666', vcount(g)), 'NODE_BORDER_COLOR'=V(g)$color), stringsAsFactors =FALSE)), 
										'edge.attributes'=list(as.data.frame(cbind(get.edgelist(g), 'EDGE_COLOR'=E(g)$color, 'EDGE_LINE_WIDTH'=E(g)$width), stringsAsFactors =FALSE)))
						colnames(graphInfo[['edge.attributes']])[1:2] <- c('NODES1', 'NODES2')
   						.convertCytoscapeWeb(graphInfo, htmlName=fileName, fileSuffix=(i-2), verbose=TRUE, destination=paste('../',fileName, '.files', sep=''), bgColor=bgColor)
						setwd(paste(fileName, '.files', sep=''))
						png(filename=tempCrossTableFileName, width=1000, height=1500, bg=bgColor)
						if (is.null(getGeneExprProfile(gAList[[i-2]]))) geneAnswersHeatmap(gAList[[i-2]], showCats=drawCats, catTerm=catTerm, geneSymbol=TRUE, nameLength=nameLength, catID=TRUE, sortBy='column', 
														colorMap=c(rgb(255-col2rgb(bgColor)[1], 255-col2rgb(bgColor)[2], 255-col2rgb(bgColor)[3], maxColorValue=255), bgColor), mapType='heatmap')
						else geneAnswersHeatmap(gAList[[i-2]], showCats=drawCats, catTerm=catTerm, geneSymbol=TRUE, nameLength=nameLength, catID=TRUE, sortBy='both', 
														colorMap=c(rgb(255-col2rgb(bgColor)[1], 255-col2rgb(bgColor)[2], 255-col2rgb(bgColor)[3], maxColorValue=255), bgColor), mapType='heatmap')
					   	#else geneAnswersConceptNet(gAList[[i-2]], centroidSize='pvalue', colorValueColumn = colorValueColumn, output='fixed', showCats=drawCats, catTerm=catTerm, geneSymbol=TRUE)
		   				dev.off()
						setwd('..')
						cat(paste('<H2 align=center><font face="courier" size="2"><A name="', paste('h-', i, sep=''), '">', tableNames[i], "</a></font></H2>", sep=''), file = outFile, sep = "\n") 
						cat('<body><center><div id="cytoscapeweb',(i-2), '" style="width:1000px;height:1000px;">Cytoscape Web will replace the contents of this div with your graph.</div></center></body>\n', file = outFile, sep="")
						cat('<center><IMG src="', paste(fileName, '.files', sep=''), '/', tempCrossTableFileName, '"></center>\n', file = outFile, sep="")
					} else {
						tempConceptFileName <- paste(tableNames[i], '_', getCategoryType(gAList[[i-2]]), '.png', sep='')
						setwd(paste(fileName, '.files', sep=''))
						png(filename=tempConceptFileName, width=1000, height=1000, bg=bgColor)
						geneAnswersConceptNet(gAList[[i-2]], centroidSize='pvalue', colorValueColumn = colorValueColumn[i-2], output='fixed', showCats=drawCats, catTerm=catTerm, geneSymbol=TRUE, bgColor=bgColor, symmetry=FALSE)
					   	#else geneAnswersConceptNet(gAList[[i-2]], centroidSize='pvalue', colorValueColumn = colorValueColumn, output='fixed', showCats=drawCats, catTerm=catTerm, geneSymbol=TRUE)
		   				dev.off()
						png(filename=tempCrossTableFileName, width=1000, height=1500, bg=bgColor)
						if (is.null(getGeneExprProfile(gAList[[i-2]]))) geneAnswersHeatmap(gAList[[i-2]], showCats=drawCats, catTerm=catTerm, geneSymbol=TRUE, nameLength=nameLength, catID=TRUE, sortBy='column', 
														colorMap=c(rgb(255-col2rgb(bgColor)[1], 255-col2rgb(bgColor)[2], 255-col2rgb(bgColor)[3], maxColorValue=255), bgColor), mapType='heatmap')
						else geneAnswersHeatmap(gAList[[i-2]], showCats=drawCats, catTerm=catTerm, geneSymbol=TRUE, nameLength=nameLength, catID=TRUE, sortBy='both', 
														colorMap=c(rgb(255-col2rgb(bgColor)[1], 255-col2rgb(bgColor)[2], 255-col2rgb(bgColor)[3], maxColorValue=255), bgColor), mapType='heatmap')
					   	#else geneAnswersConceptNet(gAList[[i-2]], centroidSize='pvalue', colorValueColumn = colorValueColumn, output='fixed', showCats=drawCats, catTerm=catTerm, geneSymbol=TRUE)
		   				dev.off()
						setwd('..') 
						cat(paste('<H2 align=center><font face="courier" size="2"><A name="', paste('h-', i, sep=''), '">', tableNames[i], "</a></font></H2>", sep=''), file = outFile, sep = "\n") 
						cat('<center><IMG src="', paste(fileName, '.files', sep=''), '/', tempConceptFileName, '"></center>\n', file = outFile, sep="")
						cat('<center><IMG src="', paste(fileName, '.files', sep=''), '/', tempCrossTableFileName, '"></center>\n', file = outFile, sep="")
						
					}
				} else {
					print('Given categories are not statistical significant!!! No category is selected!')
					cat(paste('<H2 align=center><font face="courier" size="2"><A name="', paste('h-', i, sep=''), '">', tableNames[i], "</a></font></H2>", sep=''), file = outFile, sep = "\n")
					cat('<p style="font-family:courier;text-align:center;color:red">Given categories are not statistical significant!!!!</p>', file = outFile, sep = "\n")
				}
			}
		}
		cat('<a HREF="#h-1"><p style="font-family:courier;text-align:right;font-size:11px">To Concepts-Genes Table</p></a>\n<br> <br /><br> <br />\n',file = outFile, sep = "")
	}
	
	cat('</body>', "</html>", sep = "\n", file = outFile)
	if (interactive) {
		file.remove('NetworkTransformer.jar')
		if (!keepCytoscapeFiles) {
			if (file.remove(paste(fileName, 'cytoscapeWebFiles', sep='.'))) print(paste(paste(fileName, 'cytoscapeWebFiles', sep='.'), 'is removed.'))
			else print(paste('Warning:', paste(getwd(), '', fileName, 'cytoscapeWebFiles', sep='.'), 'can not be removed.'))
		} 
	}
	close(outFile)
	print(paste(fileName, 'is generated at', getwd()))
}