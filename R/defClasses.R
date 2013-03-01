#################################################
# class def biomvRseg
#################################################
#setClass("biomvRseg",
#   representation(
#      x = "matrix",
#      segStart = "list",
#      segMean = "list",
#      group = "character",
#      family = "character"
#   ),
#   prototype = list(
#      x = matrix(),
#      segStart = list(),
#      segMean = list(),
#      group = character(),
#      family = character()
#   )
#)


### Show method
#setMethod("show", "biomvRseg", 
#	function(object) {
#		res = c(sprintf("Object is of class: 'biomvRseg'\n"),
#					sprintf("Object's data matrix is of exponential family: %s\n", object@family),
#					sprintf("Data matrix: %d x %d\n", nrow(object@x), ncol(object@x)),
#					sprintf("Column vectors belong to %d groups\n", length(unique(object@group))))
#		for(i in 1:length(unique(object@group))){
#			res<-c(res, sprintf("Column vectors group %d has %d esitmated segments\n", i, unique((sapply(object@segStart, function(x) length(x))+1)[which(object@group==unique(object@group)[i])])))
#		}
#		cat(res, "\n", sep="")
#	})
#  

### Plot method
#setMethod("plot", "biomvRseg",  
#	function(x, ...) {
#		biomvRplot(x=x@x, xSeg=x@segStart, xGrp=x@group, ...)
#	})
#	
###################################################
## class def bioMvRhsmm
###################################################	
#setClass("biomvRhsmm",
#   representation(
#      x = "matrix",
#      state= "list",
#      group = "character",
#      family = "character"
#   ),
#   prototype = list(
#      x = matrix(),
#      state= list(),
#      group = character(),
#      family = character()
#   )
#)
### Show method
#setMethod("show", "biomvRhsmm", 
#	function(object) {
#		res = c(sprintf("Object is of class: 'biomvRhsmm'\n"),
#					sprintf("Object's data matrix is of exponential family: %s\n", object@family),
#					sprintf("Data matrix: %d x %d\n", nrow(object@x), ncol(object@x)),
#					sprintf("Column vectors belong to %d groups\n", length(unique(object@group))))
#		for(i in 1:length(unique(object@group))){
#			res<-c(res, sprintf("Column vectors group %d has %d esitmated states\n", i, unique(sapply(object@state, function(v) length(unique(v)))[which(object@group==unique(object@group)[i])])))
#		}
#		cat(res, "\n", sep="")
#	})
#	
### Plot method
#setMethod("plot", "biomvRhsmm",  
#	function(x, ...) {
#		biomvRplot(x=x@x, xSeg=lapply(x@state, function(v) cumsum(runLength(Rle(v)))), xGrp=x@group, ...)
#	})	
#	

##################################################
# class def biomvRCNS using GRanges
##################################################	
setClass("biomvRCNS",
	representation(
		x = "GRanges",
		res = "GRanges",
		param = "list"
	),
	prototype = list(
		x = GRanges(),
		res = GRanges(),
		param = list()
	)
)
## Show method
setMethod("show", "biomvRCNS", 
	function(object) {
		res = c(sprintf("Object is of class: 'biomvRCNS'\n"),
					sprintf("List of parameters used in the model: %s\n", paste(names(object@param), collapse=', '))
					)
		cat(res, "\n", sep="")
		show(object@res)
	})
	
## Plot method
setMethod("plot", "biomvRCNS",  
	function(x, ...) {
		for(s in unique(values(x@res)[,'SAMPLE'])){
			for(seq in as.character(unique(seqnames(x@x)))){
				biomvRGviz(exprgr=x@x[seqnames(x@x)==seq, s], seggr=x@res[values(x@res)[,'SAMPLE']==s],...) 
			}
		}
	})	
	
	

