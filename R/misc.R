##################################################
#			segmentation functions
##################################################
#simmvRsegData<-function(nr=50, nc=10, N=c(3,4), family='norm', comVar=TRUE, segDisp=FALSE){
#	set.seed(1234)
#	## input check
#	family<-match.arg(family, c('norm', 'pois', 'nbinom'))
#	if(max(N)<nr && min(N)>1 && length(N)<=nc && length(N)>=1){
#		ng<-as.integer(rep(seq(1:length(N)), each=ceiling(nc/length(N)))[1:nc])
#	} else {
#		stop('number of changes per group not in range, or more groups than data column !')
#	}
#	
#	## generate random change position
#	segs<-lapply(N, function(x) ceiling(seq(from=2, to=nr-nr/x, length.out=x)+runif(x, 1, nr/x/2))) # generate initial change point for each grp
#	segs<-lapply(ng, function(x) segs[[x]]) # expand to a nc length list, and have each ordered
#	segs<-lapply(segs, function(x) x+floor(runif(length(x))*2)) #jitter a bit
#	rege<- lapply(segs, function(x) c(x-1, nr))
#	regs<-lapply(segs, function(x) c(1, x))
#	
#	## generate mean for each segment
#	regm<-lapply(N+1, function(x) ceiling(runif(x, 5, 10)+rep(c(-1,1), x)[1:x]*runif(x, 1, 3)))
#	regm<-lapply(ng, function(x) regm[[x]])
#	
#	## generate size factor for nbinom, alpha(dispersion) = 1/regd
#	if(family=='nbinom'){
#		if(segDisp){
#			regd<-lapply(N+1, function(x) runif(x, 50, 150))
#		} else {
#			regd<-lapply(N+1, function(x) rep(runif(1, 50, 150), x))
#		}
#		regd<-lapply(ng, function(x) regd[[x]])
#	} else {
#		regd <- vector(mode="list", length=nc)
#	}
#	
#	## generate size factor for nbinom
#	if(family=='norm'){
#		if(comVar){
#			regsd<-lapply(N+1, function(x) rep(runif(1, 1, 3), x))
#		} else {
#			regsd<-lapply(N+1, function(x) runif(x, 1, 3))	
#		}
#		regsd<-lapply(ng, function(x) regsd[[x]])
#	} else {
#		regsd <- vector(mode="list", length=nc)
#	}
#	
#	## simulate according to the family
#	x<-matrix(, nrow=nr, ncol=nc)
#	for(i in seq_len(nc)){
#		x[,i]<-unlist(lapply(seq_len(length(regm[[i]])), function(x) simUniSegData(n=rege[[i]][x]-regs[[i]][x]+1, mu=regm[[i]][x], sd=regsd[[i]][x], size=regd[[i]][x], family=family)))
#	}

#	## return simulated object
#	new("biomvRseg",
#    x = x,
#    segStart = segs,
#    segMean=regm,
#    group=ng,
#    family=family)
#}


simUniSegData<-function(n, mu, sd, size, family){
		x<-switch(family,
			norm = rnorm(n, mean=mu, sd=sd),
			pois = rpois(n, lambda=mu),
			nbinom = rnbinom(n, mu=mu, size=size))
}

mat2list<-function(x, nr=NULL){
	if(!is.null(nr)) x<-matrix(x, nrow=nr)
	L = vector(mode="list", length=nrow(x))
	for(i in 1:nrow(x)) {
		L[[i]] = cbind(estimate=x[i, 0:(i-1)]) 
	}
	return(L)
}

# check grp setting, cluster if needed, otherwise treat as one group
preClustGrp<-function(x, grp=NULL, clusterm=NULL, deepSplit=4, minClusterSize=4){
	clustmethods<-c('ward','single','complete','average','mcquitty','median','centroid')
	nc<-ncol(x)
	if (!is.null(grp)) {
		if(length(grp)!=nc) 
			stop('length of grp differs from number of column in x !!!')
		if(!is.null(clusterm)) 
			warning('clustering would not be performed, using the grp input instead !!!')	
		grp<-as.character(grp)
	}	else if(!is.null(clusterm) && nc >= 4){
		# no grp and has clusterm
		if(length(find.package('dynamicTreeCut', quiet=T))==0) {
			warning("'dynamicTreeCut' is not found, clustering result ignored, treated as one group!!!")
			grp<-rep('1', nc)
		} else{
			clusterm<-match.arg(clusterm, clustmethods)
			DM <- dist(t(x))
			HC <- hclust(DM, method=clusterm)
			require(dynamicTreeCut)
			# may provide parameter input later for this function
			grp<-cutreeDynamic(dendro = HC, distM = as.matrix(DM), deepSplit = deepSplit, pamRespectsDendro = FALSE, minClusterSize = minClusterSize, verbose=0) ## could change this minclustersize to a parameter
		}
	} else {
		# no input for grouping or cluster, thus keep all in one group
		grp<-rep('1', nc)
	}
	return(grp)
}



##################################################
#  maxgap minrun algo, utilizing RLE, mostly for transcript detection...
##################################################
maxgapminrun<-function(x, xpos=NULL, xrange=NULL, cutoff=NULL, q=0.9, minrun=5, maxgap=2, splitLen=Inf){
	# optional position vector, if not present, use index as postion
	# output values, index of active intervals, plus func parameters 
	na.rm=TRUE
	if(is.null(cutoff)){
		cutoff<- as.numeric(quantile(x, q, na.rm=na.rm))
	} else if (cutoff >= max(x, na.rm=na.rm) || cutoff<=min(x, na.rm=na.rm)) {
		stop("cutoff level is not well within range!")
	}
	
	xl<- x>=cutoff
	if(all(xl) | all(!xl)) {
		warnings("nothing to be done, all above/below cutoff !")
		intStart<-intEnd<-integer(0)
		return(list(IS=intStart, IE=intEnd, CUTOFF=cutoff, MG=maxgap, MR=minrun, SPL=splitLen))		
	}
	n<-length(x)
	
	if(!is.null(xpos) && is.null(xrange)){
		if(!is.numeric(xpos)) stop("xpos is not numeric!")
		if(n!=length(xpos)) stop("the length of xpos doesnot match the length of x!")
		if(all((xpos[-1]-xpos[-n])>maxgap) | (xpos[n]-xpos[1])<minrun)
			warnings("maxgap and minrun values are not approprate !!")
	} else if (is.null(xpos) && !is.null(xrange)){
		if(class(xrange) !='IRanges') stop("xrange is not a IRange object!")
		if(n!=length(xrange)) stop("the length of xrange doesnot match the length of x!")
		if(all(start(xrange)[-1]-end(xrange)[-n] > maxgap) | (end(xrange)[n]- start(xrange)[1]) < minrun)
			warnings("maxgap and minrun values are not approprate !!")
	}
	
	
	# this version of the maxgap minrun utilize Rle
	# every 2 consecutive runs will have different sign, a merit
	xlrle<-Rle(xl)
	i<-min(which(runValue(xlrle)))
	while(i < max(which(runValue(xlrle)))){

		rv<-runValue(xlrle)
		rl<-runLength(xlrle)	
		erl<-cumsum(rl) #idx of the end of each run
		#srl<-erl-rl+1 # idx of the start of each run, not quite used...
		nr<-length(runValue(xlrle))
		if(i == nr){
			break
		} else if(	(is.null(xpos) && is.null(xrange) && rl[i+1]<=maxgap) ||
						(!is.null(xpos) && is.null(xrange) && xpos[erl[i+1]]-xpos[erl[i]]<=maxgap) ||
						(is.null(xpos) && !is.null(xrange) && start(xrange)[erl[i+1]]-end(xrange)[erl[i]]<=maxgap) && (i+1)<nr ){
			# the next F interval is shorter than maxgap, and there is one T interval still hanging 
			# the gap distance is calculated as the distance between the end of the last True in run i to the last FALSE in the run i+1, half open
			if(	(is.null(xpos) && is.null(xrange) && sum(rl[i:min(i+2, nr)])<=splitLen) ||
				(!is.null(xpos) && is.null(xrange) && xpos[erl[min(i+2, nr)]]-xpos[erl[i]] <=splitLen) ||
				(is.null(xpos) && !is.null(xrange) && end(xrange)[erl[min(i+2, nr)]]-start(xrange)[erl[i]]<=splitLen) ){ 
				#the accumulative length of this interval including this F interval is shorter than the splitLen, merge
				runValue(xlrle)[(i+1)]<-TRUE
			} else {
				#move on and ignore this F interval, no merging
				i<-i+2
			}
		} else {
			i<-i+2
		}
	}
	
	rv<-runValue(xlrle)
	rl<-runLength(xlrle)
	erl<-cumsum(rl) 
	srl<-erl-rl+1 
	nr<-length(runValue(xlrle))
	if(!is.null(xpos)){
		ri<-which(rv & xpos[erl]-xpos[srl]>= minrun)
	} else if (!is.null(xrange)){
		ri<-which(rv & end(xrange)[erl]-start(xrange)[srl]>= minrun)
	} else {
		ri<-which(rv & rl>=minrun)
	}
	if(length(ri)>0){
		intStart<-sapply(ri, function(z) sum(rl[seq_len(z-1)])+1)
		intEnd<-sapply(ri, function(z) sum(rl[seq_len(z)]))
		#?fixme, when the first position has value above cutoff and the xpos[1] has value > then minrun, then this positon will be included as 0~xpos[1] is considered a valid run.
#		if(intEnd[1]==1 && intStart[1]==1){
#			intEnd<-intEnd[-1]
#			intStart<-intStart[-1]
#		}
		## solve neighbouring probes which are both true, but too far away
		tmp<-splitFarNeighbour(intStart=intStart, intEnd=intEnd, xpos=xpos, xrange=xrange, maxgap=maxgap, minrun=minrun)
		intStart<-tmp$IS
		intEnd<-tmp$IE
	} else {
		intStart<-intEnd<-integer(0)
	}
	return(list(IS=intStart, IE=intEnd, CUTOFF=cutoff, MG=maxgap, MR=minrun, SPL=splitLen))		
}

##################################################
# helper function to split far away neighbouring items
##################################################
splitFarNeighbour<-function(intStart=NULL, intEnd=NULL, xpos=NULL, xrange=NULL, maxgap=Inf, minrun=1){
	if(!is.null(xrange) && class(xrange)=='GRanges' || class(xrange)== 'IRanges') xpos<-NULL
	if(!is.null(xpos) || !is.null(xrange) && !is.null(intStart) && !is.null(intEnd) && length(intEnd)==length(intStart)){
		if(!is.null(xpos)){
			n<-length(xpos)
			gapStart<-which((xpos[-1]-xpos[-n])>maxgap)
		} else {
			n<-length(xrange)
			gapStart<-which((start(xrange)[-1]-end(xrange)[-n])>maxgap)	
		}
		gapEnd<-gapStart+1
		gapir<-IRanges(start=gapStart, end=gapEnd)
		gapir<-reduce(gapir)
		intN<-length(intStart)
		todel<-numeric(0)
		toadd<-numeric(0)
		for(i in 1:intN){
			splitgap<-which(intStart[i]<=gapStart & intEnd[i]>=gapEnd)
			ngap<-length(splitgap)
			if(ngap>0){
				todel<-c(todel, i)
				for(g in 1:ngap){
					if(ngap>1){
						#there are more than one
						if(g==1){
							# the 1st gap in this interval
							if(	(!is.null(xpos) && xpos[gapStart[splitgap[g]]]-xpos[intStart[i]]>=minrun) ||
								(!is.null(xrange) && end(xrange)[gapStart[splitgap[g]]]-start(xrange)[intStart[i]]>=minrun) ){
								#split first half	
								toadd<-c(toadd, intStart[i], gapStart[splitgap[g]])
							}
						}
						if (g==ngap){
							# the last gap in this interval
							if(	(!is.null(xpos) && xpos[intEnd[i]]-xpos[gapEnd[splitgap[g]]]>=minrun) ||
								(!is.null(xrange) && end(xrange)[intEnd[i]]-start(xrange)[gapEnd[splitgap[g]]]>=minrun) ){
								#split second half
								toadd<-c(toadd, gapEnd[splitgap[g]], intEnd[i])
							}
						} else{
							if(	(g+1)<ngap &&	(!is.null(xpos) && xpos[gapStart[splitgap[g+1]]]-xpos[gapEnd[splitgap[g]]]>=minrun) ||
								(!is.null(xrange) && end(xrange)[gapStart[splitgap[g+1]]]-start(xrange)[gapEnd[splitgap[g]]]>=minrun) ){
								#split second half
								toadd<-c(toadd, gapEnd[splitgap[g]], gapStart[splitgap[g+1]])
							}
						}
					} else {
						# there is only one gap within this interval
						if(	(!is.null(xpos) && xpos[gapStart[g]]-xpos[intStart[i]]>=minrun) ||
							(!is.null(xrange) && end(xrange)[gapStart[g]]-start(xrange)[intStart[i]]>=minrun) ){
							#split first half	
							toadd<-c(toadd, intStart[i], gapStart[splitgap[g]])
						}
						if( (!is.null(xpos) && xpos[intEnd[i]]-xpos[gapEnd[splitgap[g]]]>=minrun) ||
							(!is.null(xrange) && end(xrange)[intEnd[i]]-start(xrange)[gapEnd[splitgap[g]]]>=minrun) ){
							#split second half
							toadd<-c(toadd, gapEnd[splitgap[g]], intEnd[i])
						}
					}							
				}
			}
		}
		if(length(todel)>0){
			intStart<-intStart[-todel]
			intEnd<-intEnd[-todel]
		}
		if(length(toadd)>0){
			intStart<-c(intStart, toadd[seq(1, length(toadd), 2)])
			intEnd<-c(intEnd, toadd[seq(2, length(toadd), 2)])
		}	
	}
	return(list(IS=intStart, IE=intEnd))
}

gammaFit<-function(x, wt=NULL){
	if(is.null(wt)) wt <- rep(1,length(x))
	if(length(x) != length(wt)) stop("length of x and wt differ!")
	
	tmp <- cov.wt(data.frame(x),wt=wt)
    shape0 <- (tmp$center/sqrt(tmp$cov))^2
    scale<-as.numeric(tmp$center)
    shape<-optimize(function(shape) sum(dgamma(x,shape=shape, scale=scale,  log=TRUE)*wt), c(shape0-sqrt(tmp$cov), shape0+sqrt(tmp$cov)) , maximum = TRUE)[[1]]
	return(c(shape=shape, scale=scale))
}

poisFit <- function(x, wt=NULL, maxshift=1) {  	
	if(maxshift>min(x)) stop("maxshift can't be greater than the minimum of x !")
	if(is.null(wt)) wt <- rep(1/length(x),length(x))
	if(length(x) != length(wt)) stop("length of x and wt differ!")
	
	shift<-which.max(sapply(1:maxshift, function(s) dpois(x = x-s, lambda=(x-s) %*% wt,log=TRUE) %*% wt))
	lambda <- (x-shift) %*% wt 
    return(c(shift=shift, lambda=lambda))
}

nbinomFit <- function(x, wt=NULL, maxshift=1) {  	
	if(maxshift>min(x)) stop("maxshift can't be greater than the minimum of x !")
	if(is.null(wt)) wt <- rep(1,length(x))
	if(length(x) != length(wt)) stop("length of x and wt differ!")
	
	shift <- which.max(sapply(1:maxshift,  function(s) nbinomCLLDD(x, wt, s)$value))
    res<-c(nbinomCLLDD(x, wt, shift)$par, shift)
    names(res)<-c('size', 'mu', 'shift')
    return(res)
}

nbinomCLLDD<-function(x, wt=NULL, s=1){
	if(is.null(wt)) wt <- rep(1,length(x))
	if(length(x) != length(wt)) stop("length of x and wt differ!")
	m <- weighted.mean(x-s,wt)
	v <- as.numeric(cov.wt(data.frame(x-s),wt=wt)$cov)
	size <- if (v > m) m^2/(v - m) else 100
	optim(c(size,m),function(par) sum(dnbinom(x-s,size=par[1],mu=par[2],log=TRUE)*wt)  ,control=list(fnscale=-1,reltol=1e-4,maxit=1000))
}
