##################################################
# class def bioMvRseg
##################################################
setClass("biomvRseg",
   representation(
      x = "matrix",
      segStart = "list",
      segMean = "list",
      group = "integer",
      family = "character"
   ),
   prototype = list(
      x = matrix(),
      segStart = list(),
      segMean = list(),
      group = integer(0),
      family = character(0)
   )
)


## Show method
setMethod("show", "biomvRseg", 
	function(object) {
		res = c(sprintf("Object is of class: 'biomvRseg'\n"),
					sprintf("Object's data matrix is of exponential family: %s\n", object@family),
					sprintf("Data matrix: %d x %d\n", nrow(object@x), ncol(object@x)),
					sprintf("Column vectors belongs to %d groups\n", length(unique(object@group))))
		for(i in 1:length(unique(object@group))){
			res<-c(res, sprintf("Column vectors group %d has %d esitmated segments\n", i, unique((sapply(object@segStart, function(x) length(x))+1)[which(object@group==unique(object@group)[i])])))
		}
		cat(res, "\n", sep="")
	})
  

## Plot method
setMethod("plot", "biomvRseg",  
	function(x, ...) {
		biomvRplot(x=x@x, xSeg=x@segStart, xGrp=x@group, ...)
	})
	
##################################################
# class def bioMvRhsmm
##################################################	
setClass("biomvRhsmm",
   representation(
      x = "matrix",
      state= "list",
      group = "integer",
      family = "character"
   ),
   prototype = list(
      x = matrix(),
      state= list(),
      group = integer(),
      family = character(0)
   )
)
## Show method
setMethod("show", "biomvRhsmm", 
	function(object) {
		res = c(sprintf("Object is of class: 'bioMvRhsmm'\n"),
					sprintf("Object's data matrix is of exponential family: %s\n", object@family),
					sprintf("Data matrix: %d x %d\n", nrow(object@x), ncol(object@x)),
					sprintf("Column vectors belongs to %d groups\n", length(unique(object@group))))
		for(i in 1:length(unique(object@group))){
			res<-c(res, sprintf("Column vectors group %d has %d esitmated states\n", i, unique(sapply(object@state, function(v) length(unique(v)))[which(object@group==unique(object@group)[i])])))
		}
		cat(res, "\n", sep="")
	})
	
## Plot method
setMethod("plot", "biomvRhsmm",  
	function(x, ...) {
		biomvRplot(x=x@x, xSeg=lapply(x@state, function(v) cumsum(runLength(Rle(v)))), xGrp=x@group, ...)
	})	
