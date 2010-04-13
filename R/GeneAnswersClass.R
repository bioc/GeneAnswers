setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("dataframeOrNULL", c("data.frame", 'NULL'))
#---------------------------------------------------------------
# GeneAnswers: a class containing the necessary information for 
# category test of the input gene list with relative data. 
#---------------------------------------------------------------
setClass("GeneAnswers",
	representation(
		geneInput		= "data.frame",
		testType		= "character",
		pvalueT			= "numeric",
		genesInCategory	= "list",
		geneExprProfile = "dataframeOrNULL",
		annLib			= "characterOrNULL",
		categoryType	= "characterOrNULL",
		enrichmentInfo	= "data.frame" 
	),
	prototype	= list(
		geneInput		= new("data.frame"),
		testType		= "",
		pvalueT			= 1.0,
		genesInCategory	= list(),
		annLib			= "",
		categoryType	= "",
		geneExprProfile = new("data.frame"),
		enrichmentInfo  = new("data.frame")
	)
)

#---------------------------------------------------------------
# methods
#---------------------------------------------------------------
if (is.null(getGeneric("getGeneInput"))) setGeneric("getGeneInput", function(object) standardGeneric("getGeneInput"))
if (is.null(getGeneric("getTestType"))) setGeneric("getTestType", function(object) standardGeneric("getTestType")) 
if (is.null(getGeneric("getPValueT"))) setGeneric("getPValueT", function(object) standardGeneric("getPValueT")) 
if (is.null(getGeneric("getGenesInCategory"))) setGeneric("getGenesInCategory", function(object) standardGeneric("getGenesInCategory")) 
if (is.null(getGeneric("getGeneExprProfile"))) setGeneric("getGeneExprProfile", function(object) standardGeneric("getGeneExprProfile")) 
if (is.null(getGeneric("getAnnLib"))) setGeneric("getAnnLib", function(object) standardGeneric("getAnnLib")) 
if (is.null(getGeneric("getCategoryType"))) setGeneric("getCategoryType", function(object) standardGeneric("getCategoryType")) 
if (is.null(getGeneric("getEnrichmentInfo"))) setGeneric("getEnrichmentInfo", function(object) standardGeneric("getEnrichmentInfo"))    

if (is.null(getGeneric("setGeneInput"))) setGeneric("setGeneInput", function(object, geneInput) standardGeneric("setGeneInput"))
if (is.null(getGeneric("setTestType"))) setGeneric("setTestType", function(object, type=c('hyperG', 'none')) standardGeneric("setTestType")) 
if (is.null(getGeneric("setPValueT"))) setGeneric("setPValueT", function(object, pvalueT) standardGeneric("setPValueT")) 
if (is.null(getGeneric("setGenesInCategory"))) setGeneric("setGenesInCategory", function(object, genesInCategory) standardGeneric("setGenesInCategory")) 
if (is.null(getGeneric("setEnrichmentInfo"))) setGeneric("setEnrichmentInfo", function(object, enrichmentInfo) standardGeneric("setEnrichmentInfo"))
if (is.null(getGeneric("setGeneExprProfile"))) setGeneric("setGeneExprProfile", function(object, geneExprProfile) standardGeneric("setGeneExprProfile")) 
if (is.null(getGeneric("setAnnLib"))) setGeneric("setAnnLib", function(object, annLib) standardGeneric("setAnnLib")) 
if (is.null(getGeneric("setCategoryType"))) setGeneric("setCategoryType", function(object, type=c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'DOLite', 'KEGG', 'User defiend')) standardGeneric("setCategoryType")) 

if (is.null(getGeneric("summary"))) setGeneric("summary", function(object) standardGeneric("summary"))
if (is.null(getGeneric("show"))) setGeneric("show", function(object) standardGeneric("show"))

setMethod("getGeneInput",signature(object="GeneAnswers"), function(object) object@geneInput)

setMethod("getTestType",signature(object="GeneAnswers"), function(object) object@testType)

setMethod("getPValueT",signature(object="GeneAnswers"), function(object) object@pvalueT)

setMethod("getGenesInCategory",signature(object="GeneAnswers"), function(object) object@genesInCategory)

setMethod("getGeneExprProfile",signature(object="GeneAnswers"), function(object) object@geneExprProfile)

setMethod("getAnnLib",signature(object="GeneAnswers"), function(object) object@annLib)

setMethod("getCategoryType",signature(object="GeneAnswers"), function(object) object@categoryType)

setMethod("getEnrichmentInfo",signature(object="GeneAnswers"), function(object) object@enrichmentInfo)

setMethod("setGeneInput", signature(object="GeneAnswers"), function(object, geneInput)
{
	if (is.vector(geneInput) | is.data.frame(geneInput) | is.matrix(geneInput)) {
		object@geneInput <- geneInput
		return(object)
	} else {
		stop('Input is not a vector or data frame or matrix! No value is set')
	}
})

setMethod("setTestType", signature(object="GeneAnswers"), function(object, type=c('hyperG', 'none'))
{
	type <- match.arg(type)
	if ((type == 'none') & !is.null(object@enrichmentInfo)) stop('Enrichment has existed! No value is set')
	object@testType <- type
	return(object)
})

setMethod("setPValueT", signature(object="GeneAnswers"), function(object, pvalueT)
{
	if (is.numeric(pvalueT) & (pvalueT <= 1) & (pvalueT > 0)) {
		object@pvalueT <- pvalueT
		return(object)
	} else {
		stop('Input is not a valid p value threshold! No value is set')
	}
})

setMethod("setGenesInCategory", signature(object="GeneAnswers"), function(object, genesInCategory)
{
	if (is.list(genesInCategory)) {
		object@genesInCategory <- genesInCategory
		return(object)
	} else {
		stop('Input is not a list! No value is set')
	}
})

setMethod("setEnrichmentInfo", signature(object="GeneAnswers"), function(object, enrichmentInfo)
{
	if (is.data.frame(enrichmentInfo) | is.matrix(enrichmentInfo)) {
		object@genrichmentInfo <- as.data.frame(enrichmentInfo)
		return(object)
	} else {
		stop('Input is not a list! No value is set')
	}
})


setMethod("setGeneExprProfile", signature(object="GeneAnswers"), function(object, geneExprProfile)
{
	if (is.data.frame(geneExprProfile) | is.matrix(geneExprProfile)) {
		object@ggeneExprProfile <- geneExprProfile
		return(object)
	} else {
		stop('Input is not a data frame or matrix! No value is set')
	}	
})

setMethod("setAnnLib", signature(object="GeneAnswers"), function(object, annLib)
{
	if (is.character(annLib) | is.null(annLib)) { 
		object@annLib <- annLib
		if (is.null(annLib)) object@categoryType <- 'User defiend'
		return(object)
	} else {
		stop('Input is not a valid annotation libraray! No value is set')
	}
})

setMethod("setCategoryType", signature(object="GeneAnswers"), function(object, type=c('GO', 'GO.BP', 'GO.CC', 'GO.MF', 'DOLite', 'KEGG', 'User defiend'))
{
	type <- match.arg(type)
	object@categoryType <- type
	if (type == 'User defiend') object@annLib <- NULL
	return(object)
})

setMethod("summary",signature(object="GeneAnswers"), function(object)
{
	cat(paste('This GeneAnswers instance was build from', object@categoryType, 'based on', object@testType, 'test.\n'))
	cat(paste('Statistical information of', dim(object@enrichmentInfo)[1], 'categories with p value less than', object@pvalueT, 'are reported. Other categories are considered as nonsignificant.\n'))
	cat(paste('There are', length(object@genesInCategory), 'categories related to the given', dim(object@geneInput)[1], 'genes\n'))
	cat('\n')
	show(object)
})

setMethod("show",signature(object="GeneAnswers"), function(object)
{
	cat('Summary of GeneAnswers instance information:\n')
	cat('\nSlot: geneInput\n')
	print(head(object@geneInput))
	if (dim(object@geneInput)[1] > 6) cat('......\n')
	cat('\nSlot: testType\n')
	print(object@testType)
	cat('\nSlot: pvalueT\n')
	print(object@pvalueT)
	cat('\nSlot: genesInCategory\n')
	print(head(object@genesInCategory))
	if (length(object@genesInCategory) > 6) cat('......\n')
	cat('\nSlot: geneExprProfile\n')
	print(head(object@geneExprProfile))
	if (!is.null(object@geneExprProfile)) {
		if (dim(object@geneExprProfile)[1] > 6) cat('......\n')
	}
	cat('\nSlot: annLib\n')
	print(object@annLib)
	cat('\nSlot: categoryType\n')
	print(object@categoryType)
	cat('\nSlot: enrichmentInfo\n')
	print(head(object@enrichmentInfo))
	if (dim(object@enrichmentInfo)[1] > 6) cat('......\n')
})