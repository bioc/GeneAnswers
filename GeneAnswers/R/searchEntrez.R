`searchEntrez` <-
function(tagList, species='human') {
	return(lapply(tagList, .searchEntrezTag, species=species))
}

