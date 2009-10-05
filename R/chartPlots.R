`chartPlots` <-
function(x, chartType=c('pieChart', 'barPlot', 'all'), specifiedCols=c('genes in Category'), top=5, newWindow=TRUE, ...) {
	chartType <- match.arg(chartType)
	if (is.matrix(x) | is.data.frame(x)) {
		if (!(FALSE %in% (specifiedCols %in% colnames(x))) & (length(specifiedCols) == 1)) {
			if ((chartType == 'pieChart') | (chartType == 'all')) {
				if (newWindow) quartz()
#				if (toupper(Sys.info()[1]) == 'WINDOWS') windows()
#				else quartz()
				pie(as.numeric(x[1:top, specifiedCols[1]]), labels=rownames(x)[1:top], col = rainbow(top), radius = 0.9, main=paste('Top ', top, ' Categories Distribution based on ', specifiedCols, sep=''), ...)
			}
			if ((chartType == 'barPlot') | (chartType == 'all')) {
				if (newWindow | chartType == 'all') quartz()
#				if (toupper(Sys.info()[1]) == 'WINDOWS') windows()
#				else quartz()
				barplot(as.numeric(x[1:top, specifiedCols[1]]), names.arg =rownames(x)[1:top], space=0.6, col=rainbow(5), main=paste('Top ', top, ' Categories Distribution based on ', specifiedCols, sep=''), ...)
			}
		} else {
			stop('One or two given names are not column names in the original matrix!')
		}
	} else {
		stop('The input data format is not a matrix or dataframe!')
	}
}

