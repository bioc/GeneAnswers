`getListGIF` <- function(glist=glist, output=NULL, background=c("white", "black"), type=c("genelevel", "riflevel", "goterm")) {
	background <- match.arg(background)
	type <- match.arg(type)
	if (is.null(glist)) stop('Please provide a list of gene symbols')
	if (!is.vector(glist)){
		if(is.list(glist)){
			glist <- unlist(glist)
		} else if(is.data.frame(glist)) {
			glist <- glist[,,drop = TRUE]
		}
	}
	
	if(length(glist)==1){
		querylist <- glist
	} else {
		## concatenate the gene list with coma as the delimiter
		querylist <- paste(glist, collapse = ',')
	}
	
	if(!(background %in% c('white', 'black'))) {
		message('Background must be either "white" or "black"!')
		background <- 'white'
	}
	
	if(!(type %in% c('genelevel', 'riflevel', 'goterm'))) {
		stop('Type must be either "genelevel", "riflevel" or "goterm"!')
	}
	
	## create url to submit to ListGIF server  
	submit.url <- paste('http://listgif.nubic.northwestern.edu/webservice.php?BG=', background, '&type=', type, '&list=', querylist, sep='') 
	html <- htmlTreeParse(submit.url, useInternalNodes=TRUE)
	url.png <- xpathSApply(html, "//*/div", xmlValue)
	url.png <- gsub("\\n|\\t", "", url.png)
	if(!is.null(output)){
		#download.file(url.png, destfile=filename)
		result <- try(download(url.png, destfile=output, mode = 'wb', quiet=TRUE))
	} else{ 
		png.default <- unlist(strsplit(url.png, '/'))
		png.default <- png.default[length(png.default)]
		result <- try(download(url.png, destfile=png.default, mode = 'w', quiet=TRUE))
    }
}
