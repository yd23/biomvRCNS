##################################################
# retain non-parametric, possion, gamma, and nbinom for the sojourn dist, these could be directly fork from mhsmm
##################################################

##################################################
# re-implement HSMM, one more slot to handle distance array
##################################################
biomvRhsmm<-function(x, maxk=NULL, maxbp=NULL, J=3, xPos=NULL, xRange=NULL, usePos='start', xAnno=NULL, soj.type='gamma', emis.type='norm', q.alpha=0.05, r.var=0.75, iterative=TRUE, maxit=1000, maxgap=Inf, tol=1e-04, grp=NULL, clusterm=NULL, na.rm=TRUE){
	## input checking
	# x, matrix/range like
	# xPos, x feature information if x is not a grange object; so for count data, x should be a GRange , for other continous data, a matrix /  Grange
	# maxk, maximum possible stays of a states; maxbp, maximum possible duration of stay in bp, give xPos has the coordinates.
	# J, number of states
	# xAnno, dataframe/range like object contains the same number of features for x data rows, then consider to lock sojourn dist
	# soj.type, sojourn dist type
	# family, whether a array like normal/t, or sequencing like poisson/nbinom
	# iterative, default T for count, and F for real
	# maxit, max itration 
	# ints, prior information, need more, see core.r	#fixme
	# lock.transition / lock.d, lock transition and sojourn #fixme
	# est.method=c('viterbi', 'smooth')
	
	#  tol=1e-4
	#  ksmooth.thresh = 1e-20 #this is a threshold for which d(u) values to use - if we throw too many weights in the default density() seems to work quite poorly
	#  shiftthresh = 1e-20 #threshold for effective "0" when considering d(u)	
	
	if (!is.numeric(x) &&  !is.matrix(x) && class(x)!='GRanges') 
        stop("'x' must be a numeric vector or matrix or a GRanges object.")
    if(class(x)=='GRanges') {
    	xid<-names(values(x))
    	xRange<-x
    	mcols(xRange)<-NULL
    	x<-as.matrix(values(x))
    } else if(length(dim(x))==2){
		xid<-colnames(x)
		x<-as.matrix(x)
	} else {
		warning('no dim attributes, coercing x to a matrix with 1 column !!!')
		x <- matrix(as.numeric(x), ncol=1)
		xid<-paste('S', seq_len(nc), sep='')
		colnames(x)<-xid
	}
	nr<-nrow(x) 
	nc<-ncol(x)
	
	## some checking on xpos and xrange, xrange exist then xpos drived from xrange,
	if(!is.null(xRange) && (class(xRange)=='GRanges' || class(xRange)=='IRanges') && !is.null(usePos) && length(xRange)==nr && usePos %in% c('start', 'end', 'mid')){
		if(usePos=='start'){
			xPos<-start(xRange)
		} else if(usePos=='end'){
			xPos<-end(xRange)
		} else {
			xPos<-(start(xRange)+end(xRange))/2
		}
	} else {
		# no valid xRange, set it to null
		warning('no valid xRange and usePos found, check if you have specified xRange / usePos.')
		xRange<- NULL
	} 
	if (is.null(xPos) || !is.numeric(xPos) || length(xPos)!=nr){
		warnings("No valid positional information found. Re-check if you have specified any xPos / xRange.")
		xPos<-NULL
	}
	if (!is.null(maxbp) && (!is.numeric(maxbp) || (length(maxbp) != 1) || (maxbp <= 1) ||  ( !is.null(xPos) && maxbp > max(xPos,na.rm=na.rm)-min(xPos, na.rm=na.rm)))) 
	 	 stop(sprintf("'maxk' must be a single integer between 2 and the maximum length of the region if xPos is avaliable."))	
	
	# check grp setting, cluster if needed, otherwise treat as one group
	if(!is.null(grp)) grp<-as.character(grp)
	grp<-preClustGrp(x, grp=grp, clusterm=clusterm)
	
	# initial sojourn setup unifiy parameter input / density input,  using extra distance, non-integer value can give a dtype value
	if(!is.null(xAnno) && !is.null(soj.type) && soj.type %in% c('gamma') && class(xAnno) %in% c('TranscriptDb', 'GRanges', 'GRangesList')){
		#	this is only used when the xAnno object contains appropriate annotation infromation which could be used as prior for the sojourn dist in the new HSMM model
		# if xAnno is also present, then J will be estimated from xAnno, and pop a warning, ## this only make sense if difference exist in the distribution of sojourn of states.		
		soj<-sojournAnno(xAnno, soj.type=soj.type)
		J<-soj$J
		cat('Estimated state number from xAnno: J = ', J, '\n', sep='')	
		#now, if there is xAnno, then J, maxbp could be infered, and if xPos exists, then maxk could also be infered.
		if(is.null(maxbp)){
			# estimating a reasonable number for maxbp
#			fixme, this is related to the soj type, now only have gamma avaliable
			maxbp<-ceiling(median(sapply(1:J, function(j) soj$shape[j]*soj$scale[j])))	
		}
		soj<-append(soj, maxbp=maxbp) 
	} else if(is.numeric(J) && J>1){
		cat('xAnno is not present or not supported, try to use maxbp/maxk in the uniform prior for the sojourn distribution!\n') 
		# J is ok
		if(is.null(xPos)){
			# no position as well, in turn means no xRange nor multiple seq, 
			# init if J and maxk are ok
			if (!is.null(maxk) && is.numeric(maxk) && (length(maxk) == 1) && (maxk > 1) &&  (maxk < nr)) {
				warning('maxbp and xPos are not present or not valid, using maxk for the sojourn distribution !!!')
				soj<-list(d = unifMJ(maxk, J), type = soj.type, J=J, maxk=maxk)
			} else {
				stop(sprintf("'maxk' must be a single integer between 2 and the number of rows of 'x': %d.", nr))
			}
		} else if(!is.null(maxbp) && maxbp > 1){
			# has position and good maxbp, will init it latter, maxk will be estimated there
			soj<-list(J=J, maxbp=maxbp, type = soj.type)
		} else {
			stop(sprintf("Both maxk and maxbp are not avaliable!"))
		}
	} else {
		# no good J
		stop("J must be specified or estimated from xAnno !!")
	}
	# so far, soj is a list object, depending on which case
	# case1, soj param J, maxbp from xAnno
	# case2, no pos, has input maxk and J and d, ready for initSojDd
	# case3, has input pos and maxbp

	

	# check mv vs iterative
	if(emis.type=='mvnorm' && iterative==TRUE){
		emis.type<-'norm'
	} else if (emis.type=='mvnorm' ){
		iterative<-FALSE
	}
	
	## build xRange if not a GRanges for the returning object
	if(is.null(xRange) || class(xRange) != 'GRanges'){
		if(!is.null(xRange) && class(xRange) == 'IRanges'){
			xRange<-GRanges(seqnames='biomvRCNS', xRange)	
		} else 	if(!is.null(xPos)){
			xRange<-GRanges(seqnames='biomvRCNS', IRanges(start=xPos, width=1))	
		} else {
			xRange<-GRanges(seqnames='biomvRCNS', IRanges(start=seq_len(nr), width=1))	
		}
	}
	# get seqnames status	
	seqs<-unique(as.character(seqnames(xRange)))
	
	
	## initialize the output vectors
	state <- matrix(NA, nrow=nr, ncol=nc)
	res<-GRanges(); #seqlevels(res)<-seqlevels(xRange)
	for(g in unique(grp)){
		cat(sprintf("step 1 building HSMM for group %s\n", g))
		gi<-grp==g
		if(length(seqs)>1){
			# we have more than one seq to batch
			# need to 'c' res into resl[[c]]
			for(s in seq_len(length(seqs))){
				r<-as.character(seqnames(xRange)) == seqs[s]
				# prep soj for the c loop, since there are multiple seq, which also means there must be xpos and maxbp
				ssoj<-append(soj, initDposV(xPos[r], maxbp))
				if(is.null(ssoj$fttypes)){
					ssoj<-append(ssoj, list(d=unifMJ(ssoj$maxk*sum(r), J)))
				}
				ssoj <- initSojDd(ssoj)
				
				if(iterative){
					for(c in which(gi)){
						cat(sprintf("step 1 building HSMM for seq %s in column %s.\n", seqs[s], c))
						runout<-biomvRhsmmRun(x[r,c], xid[c], xRange[r], ssoj, emis.type, q.alpha, r.var, maxit)	
						res<-c(res, runout$res)
						state[r, c]<-runout$yhat
					}
				} else {
					runout<-biomvRhsmmRun(x[r,gi], xid[gi], xRange[r], ssoj, emis.type, q.alpha, r.var, maxit)	
					res<-c(res, runout$res)
					state[r, gi]<-runout$yhat
				}
			}
		} else {
			# only one seq present
			# prep soj for all
			if(is.null(soj$d)){
				# not case 2
				ssoj<-append(soj, initDposV(xPos, maxbp))
				if(is.null(soj$fttypes)){
					ssoj<-append(ssoj, list(d=unifMJ(ssoj$maxk*nr, J)))
				}
			}
			ssoj <- initSojDd(ssoj)
			if(iterative){
				for(c in which(gi)){
					cat(sprintf("step 1 building HSMM for column %s\n", c))
					runout<-biomvRhsmmRun(x[,c], xid[c], xRange, ssoj,emis.type, q.alpha, r.var, maxit)
					res<-c(res, runout$res)
					state[,c]<-runout$yhat
				}
			} else {
				runout<-biomvRhsmmRun(x[,gi], xid[gi], xRange, ssoj, emis.type, q.alpha, r.var, maxit)
				res<-c(res, runout$res)
				state[,gi]<-lapply(seq_len(sum(gi)), function(c) runout$yhat)	
			}			
		}	
	} # end for g
	
	
	# setup input data and state to xRange for returning
	colnames(state)<-paste('state.',xid, sep='')
	values(xRange)<-DataFrame(x, state, row.names = NULL)
	new("biomvRCNS",  
		x = xRange, res = res,
		param=list(J=J, maxk=maxk, maxbp=maxbp, maxgap=maxgap, soj.type=soj.type, emis.type=emis.type, q.alpha=q.alpha, r.var=r.var, iterative=iterative, maxit=maxit, tol=tol, group=grp, clusterm=clusterm, na.rm=na.rm)
	)
}



biomvRhsmmRun<-function(x, xid='sampleid', xRange, soj, emis.type='norm', q.alpha=0.05, r.var=0.75, maxit=1, na.rm=TRUE){
	# now x should be a one column matrix
	if(is.null(dim(x))) x<-matrix(x); colnames(x)<-xid
	nr<-nrow(x)
	J<-soj$J
	maxk<-soj$maxk
	
	
	# create default uniform initial probablity
	init<-rep(1/J, J) # start with uniform
	# create default uniform transition probablity
	trans <- matrix(1/(J-1), nrow = J, ncol=J)
	diag(trans)<-0
	
	# initialize emission parameters, either from user input or raw data
	emis<-list(type=emis.type)
	if(emis$type == 'norm' || emis$type== 'mvnorm') {
		emis$mu <- estEmisMu(x, J, q.alpha=q.alpha)
		emis$var <- estEmisVar(x, J, r.var=r.var)
	} else if (emis$type == 'pois'){
		emis$lambda  <- estEmisMu(x, J, q.alpha=q.alpha)
	}
	
	#estimation of most likely state sequence
	# switch est.method   viterbi , .C /  smooth
	#define likelihood
	ll <- rep(NA,maxit)
	# start MM iteration
	for(it in 1:maxit) {
		# reestimationg of emmision   
		emis<-initEmis(emis=emis, x=x)
		B  = .C("backward", a=as.double(trans), pi=as.double(init), b=as.double(emis$p), d=as.double(soj$d), D=as.double(soj$D),
				  maxk=as.integer(maxk), DL=as.integer(nrow(soj$d)), T=as.integer(nr), J=as.integer(J), 
				  eta = double(nrow(soj$d)*J), L=double(nr*J), N=double(nr), ahat=double(J*J), pihat=double(J),
				  F=double(nr*J), G=double(nr*J), L1 = double(nr*J), si=double(nr*J))#, PACKAGE='bioCNS')

		#check gamma and eta
#		if(any(is.nan(B$L))) {
#		  warnings("NaNs detected in gamma.  Exiting...")
#		  return(B)
#		}
		if(any(B$L<0)) B$L = zapsmall(B$L)      
#					if(any(B$eta<0)) B$eta = zapsmall(B$eta)      
#					if(any(B$N<0))  B$N = zapsmall(B$N)		

		#update initial prob PI, transition >=0 check
		init<-B$pihat
		init[init<0]<-0
		trans <- matrix(B$ahat,ncol=J)
		trans[trans<0] <- 0
	
		#update emision according to the new estimated distribution paramenters using B$gamma sample in  (mstep(x,matrix(B$gamma,ncol=J))))	
		emis<-initEmis(emis=emis, x=x, B=B)
		# update sojourn dD, using eta(nonparametric), eta+d (ksmoothed-nonparametric), eta+d+shift(poisson), eta(nbinom), eta(gamma)
		soj<-initSojDd(soj=soj, B=B)

		# loglikelihood for this it, using B$N
		ll[it]<-sum(log(B$N))
		if( it>1 && abs(ll[it]-ll[it-1]) < tol) {
			break()	
		}
	}	 # end for maxit
	## assign states and split if necessary.
	if(is.null(soj$fttypes)){
		yhat<-as.character(apply(matrix(B$L,ncol=J),1,which.max))
	} else {
		yhat<-soj$fttypes[apply(matrix(B$L,ncol=J),1,which.max)]
	}

	# setup this new res gr
	Ilist<-lapply(unique(yhat), function(j) do.call(cbind, splitFarNeighbouryhat(yhat, xRange=ranges(xRange), maxgap=Inf, state=j)))
	names(Ilist)<-unique(yhat)
	res<- do.call('c', 
					lapply(unique(yhat), function(j)
						GRanges(seqnames=as.character(seqnames(xRange)[1]), 
							IRanges(start=rep(start(xRange)[Ilist[[j]][,'IS']], length(xid)), end=rep(end(xRange)[Ilist[[j]][,'IE']], length(xid))), 
							SAMPLE=rep(xid, each=nrow(Ilist[[j]])), 
							STATE=rep(as.character(j), nrow(Ilist[[j]])*length(xid)), 
							MEAN=as.numeric(sapply(xid, function(s) apply(Ilist[[j]], 1, function(r) apply(as.matrix(x[r[1]:r[2],s]), 2, mean, na.rm=na.rm))))
						)
					)
			)
	return(list(yhat=yhat, res=res))
}




##################################################
#			HSMM helper functions
##################################################
##################################################
# create lists of uniform prior
##################################################
unifMJ<-function(M,J, ints=NULL){
	# not finished create unif with supplied ints
	if(is.null(ints)){
		ret<-do.call(cbind, lapply(1:J, function(j) dunif(1:M, 1, M)))		
	} else if(min(ints)<1 | max(ints)> M) {
		stop('All supplied intervals should be in the range of [1, M]')
	} else {
		# ints, a dataframe marks the starts and ends for each interval
		# fixme, seems need some more thinking, currently, for all states the unif prior has the same interval setup.
		ret<-do.call(cbind, lapply(1:J, function(j) do.call(c, lapply(1:nrow(ints), function(i) dunif(1:M, ints[i,1], ints[i,2])))))		
	}
	ret
}



##################################################
# simulate annotation for simluated dataset, into 3 group, monotonicly accending x value
##################################################
simuAnno<-function(obj){
	if(class(obj) != 'biomvRseg') stop("invalid class for input object, must be 'biomvRseg' !")
	nr<-nrow(obj@x)
	# with multiple groups, then the annotation simulation should be done for each group
	# thus will generate ngroup 
	gtag<-unique(obj@group)
	ng<-length(gtag)
	if(ng>1){
		# more than one group
		xAnno<-GRangesList()
		ptype<-matrix(0, nrow=nr, ncol=ng)
		colnames(ptype)<-gtag
		for(g in 1:ng){
			s<-match(gtag[g], obj@group)
			segType<-rep('B', length(obj@segMean[[s]]))
			segType[which(obj@segMean[[s]]<quantile(obj@segMean[[s]], 0.25))]<-'A'
			segType[which(obj@segMean[[s]]>quantile(obj@segMean[[s]], 0.7))]<-'C'
			pr<-c(1, obj@segStart[[s]], nr+1)
			ptype[,g]<-unlist(sapply(2:length(pr), function(i) rep(segType[i-1], pr[i]-pr[i-1])))
			prle<-Rle(ptype[,g])
			xAnno<-c(xAnno, GRangesList(GRanges(seqnames = Rle(c("chr99"), ), ranges = IRanges(start=c(1, (cumsum(runLength(prle))+1)[1:(length(runLength(prle))-1)]), end=cumsum(runLength(prle))) , type= runValue(prle))))
		}
		
	} else {
		#single group
		segType<-rep('B', length(obj@segMean[[1]]))
		segType[which(obj@segMean[[1]]<quantile(obj@segMean[[1]], 0.25))]<-'A'
		segType[which(obj@segMean[[1]]>quantile(obj@segMean[[1]], 0.7))]<-'C'
		pr<-c(1, obj@segStart[[1]], nr+1)
		ptype<-unlist(sapply(2:length(pr), function(i) rep(segType[i-1], pr[i]-pr[i-1])))
		prle<-Rle(ptype)
		xAnno<-GRanges(seqnames = Rle(c("chr99"), ), ranges = IRanges(start=c(1, (cumsum(runLength(prle))+1)[1:(length(runLength(prle))-1)]), end=cumsum(runLength(prle))) , type= runValue(Rle(ptype)))
	}
	# when corecing a multicolumn data.fram / matrix with chars to the values() of granges, careful
	xMeta<-GRanges(seqnames = Rle(c("chr99"), nr), ranges = IRanges(start=1:nr, width=1) , type=data.frame(X1=ptype, stringsAsFactors=F))
	return(list(xMeta=xMeta, xAnno=xAnno))
}

gammafit <- function(x,wt=NULL) {
	tol = 1e-08
	if(is.null(wt)) wt = rep(1,length(x))

      tmp = cov.wt(data.frame(x),wt=wt)
      xhat = tmp$center
      xs = sqrt(tmp$cov)
      s = log(xhat) - mean(weighted.mean(log(x),wt))    
      aold = (xhat/xs)^2
      a = Inf
      while(abs(a-aold)>tol) {
        a = aold - (log(aold) - digamma(aold) - s)/((1/aold) - trigamma(aold))        
        aold=a
      }
      return(list(shape=a,scale=xhat))
}
##################################################
# estimate sojourn distribution from annotation using gamma
##################################################
sojournAnno<-function(xAnno, soj.type= 'gamma', pbdist=NULL){ 
	# xAnno has to be a Grange / rangedata obj
	# check if xAnno class, txdb or dataframe or rangedata ...
	# there is also the possiblity of proposing an emperical number for the states.
	# must ensure there are at least 2 for each state ? todo
		
	if(class(xAnno) == 'TranscriptDb') {
		J<-3
		fttypes<-c('intergenic', 'intron', 'exon')
		#3 feature type, exon, intron, intergenic
		transc <- transcripts(xAnno) # this give you all cds ranges ungroupped
		intergenic<-gaps(transc)
		exon <- exons(xAnno) # this give you all exon ranges ungroupped
		intron<- unlist(intronsByTranscript(xAnno))
		ftdist<-list(intergenic=width(intergenic), intron=width(intron), exon=width(exon))
	} else if(class(xAnno)=='GRanges'){
		# then the first column of the elementMetadata must be a character vector marking the type of features.
		# need a way to sort 
		fts<-values(xAnno)[,1]
		fttypes<- unique(fts)
		J<-length(fttypes)
		ftdist<-lapply(1:J, function(j) width(xAnno[fts==fttypes[j]]))
	} else if (class(xAnno)=='GRangesList'){
		ng<-length(xAnno)
		## the assumptions are number of ft types could be different from differnet list entry, group wise analysis, which of coz could be wraped using foreach(ng) single group approach
		fts<-lapply(1:ng, function(g) values(xAnno[[g]])[,1])
		fttypesL<- lapply(fts, function(x) unique(x[order(x)]))
		J<-length(unique(unlist(fts)))
		ftdist<-lapply(1:ng, function(g) lapply(1:length(fttypesL[[g]]), function(j) width(xAnno[[g]][fts[[g]]==fttypesL[[g]][j]])))
		fttypes<-unique(unlist(fttypesL))
		fttypes<-fttypes[order(fttypes)]
		ftdist<-lapply(fttypes, function(t) unlist(lapply(1:ng, function(g) ftdist[[g]][[match(t, fttypesL[[g]])]])))
	}
	
	# should switch here between different soj.type	
	if(soj.type=='gamma'){
		shape<-numeric()
		scale<-numeric()
		for(i in 1:J){
			# fit  for all feature type
			# should use this func, switch here between soj dist.
			gampar<-gammafit(ftdist[[i]])
			if(! is.null(pbdist)){
				# if distance between points are even
				gampar$scale <- pbdist / gampar$shape
			}
			shape<-c(shape, gampar$shape)
			scale<-c(scale, gampar$scale)
		}
		soj<-list(shape=shape, scale=scale, type = soj.type, fttypes=fttypes, J=J)
	}

	#return soj object
	return(soj)
}

initDposV<-function(xpos, maxbp){
	# for each position, find the maxk
	nr<-length(xpos)
	maxbpidx<-sapply(1:nr, function(i) max(which(xpos[i]+maxbp >= xpos)))
	# find the maxk idx
	maxk<-max(maxbpidx - seq_len(nr))+1
	# initialize the maxk position list for each position
	## option 2, a TM * J matrix
	dposV<-c(sapply(1:nr, function(t) xpos[t:(t+maxk-1)]-xpos[t]))+1 # dposV[(t-1)*maxk+u]
	#dposV[1]<-.Machine$double.eps # this will not ensure a 1 for the first position after the dgamma call.
	# sub > maxbp and NA
	dposV[which(dposV>maxbp)]<-NA
	dposV[is.na(dposV)]<-.Machine$double.xmax # for those positions exceed maxbp, this works fine #fixme, now there is another problem in b/f algo, the xmax will make the reestimation of soj$d difficult...
	return(list(dposV=dposV, maxk=maxk))
}




initSojDd <- function(soj, B=NULL) {
	# take the initial soj dist parameter / or sample of density
	if(! soj$type %in% c('nparam', 'gamma', 'poisson', 'nbinom')) stop("invalid sojourn type found in soj$type !")
	
	# parameter initialisation
	J<-soj$J
	maxk<-soj$maxk
	if(is.null(soj$dposV)){
		dposV<-1:maxk
	} else {
		## here d for all positions should be aggregated within each J
		dposV<-soj$dposV
	}
	idx<-dposV != .Machine$double.xmax
	nb <- length(dposV)/maxk
	if(soj$type == "gamma") {
		if(!is.null(B)){
			# then this is a update run
			soj$d <- matrix(B$eta+.Machine$double.eps,ncol=J)
			soj$shape <- soj$scale <- numeric(J)
			for(j in 1:J) {           
				param <- gammafit(dposV[idx],wt=soj$d[idx,j])
				soj$shape[j] <- param$shape
				soj$scale[j] <- param$scale     	  
			}
		}	          
		if(!is.null(soj$shape) && !is.null(soj$scale) && length(soj$shape)==J && length(soj$scale)==J) {
			# for update with para estimated from B, or initial sojourn using param
			soj$d<-sapply(1:J, function(j) dgamma(dposV,shape=soj$shape[j],scale=soj$scale[j]))
		} # else assume soj$d exist.
	}
	# add D slot
	soj$D <- sapply(1:J, function(j) sapply(1:nb, function(t) rev(cumsum(rev(soj$d[((t-1)*maxk+1):(t*maxk),j])))))
	# return soj with the initial d slot
	return(soj)
}

##################################################
# to estimate segment wise mean vector/list
##################################################
estEmisMu<- function(x, J, q.alpha=0.05, na.rm=TRUE){
	nc<-ncol(x)
	if(J >1){
		if(is.null(nc) || nc == 1){
			# univariate
			ret<-as.numeric(quantile(x, seq(from=q.alpha, to=(1-q.alpha), length.out=J), na.rm=na.rm))
		} else {
			# multiple
			muv<-t(sapply(1:nc, function(i) as.numeric(quantile(x[,i], 	seq(from=q.alpha, to=(1-q.alpha), length.out=J), na.rm=na.rm))))
			ret<-lapply(apply(muv, 2, list),unlist)
		}
	} else {
		ret<-mean(x, na.rm=na.rm)
	}
	return(ret)
}
##################################################
# to estimate segment wise variance vector / covariance matrix list
##################################################
estEmisVar<-function(x, J=3, na.rm=TRUE, r.var=0.75){
	# q.var is the espected ration of variance for state 1 and J versus any intermediate states
	# a value larger than 1 tend to givce more extreme states;  a value smaller than 1 will decrease the probablity of having extreme state, pushing it to the center.
	nc<-ncol(x)
	f.var<-rep(ifelse(r.var>=1, 1, r.var), J)
	if(J%%2 == 1) {
		f.var[(J+1)/2]<-ifelse(r.var>=1, 1/r.var, 1)
	} else if( J>2 ){
		f.var[(J/2):(J/2+1)]<-ifelse(r.var>=1, 1/r.var, 1)
	}	
	if(J >1){
		if(is.null(nc) || nc == 1){
			# univariate
			ret<-var(x, na.rm=na.rm)*f.var
		} else {
			# multiple
			if(na.rm) na.rm<-'complete.obs'
			ret<-lapply(1:J, function(j) cov(x, use=na.rm)*f.var[j])
		}
	} else {
		ret<-var(as.numeric(x), na.rm=na.rm)
	}
	return(ret)
}
##################################################
# initialize and update emission probablity
##################################################
initEmis<-function(emis, x, B=NULL){
	if(is.null(B)){
		# then this is for p initialization
		if(emis$type == 'mvnorm') {
			J<-length(emis$mu)
			emis$p <- sapply(1:J, function(j) dmvnorm(x, mean = emis$mu[[j]],  sigma = emis$var[[j]])) # here sigma requires covariance mat
		} else if (emis$type == 'pois'){
			J<-length(emis$lambda)
			emis$p <-sapply(1:J, function(j) dpois(x, lambda=emis$lambda[j]))
		} else if (emis$type == 'norm'){
			J<-length(emis$mu)
			emis$p <- sapply(1:J, function(j) dnorm(x,mean=emis$mu[j], sd=sqrt(emis$var[j]))) # here is sd
		}
	} else {
		# then this is for the re-restimation of emis param
		J<-B$J
		if(emis$type == 'mvnorm') {
			isa <-  !apply(is.na(x),1,any) # Find rows with NA's (cov.wt does not like them)
			tmp <- apply(matrix(B$L,ncol=J), 2, function(cv) cov.wt(x[isa, ], cv[isa])[c('cov', 'center')]) # x is already a matrix 
			emis$mu <- sapply(tmp, function(l) l['cneter'])
			emis$cov <- sapply(tmp, function(l) l['cov'])
		} else if (emis$type == 'pois'){
			emis$lambda = sapply(1:J, function(j) weighted.mean(x, B$L[((j-1)*nr+1):(j*nr)]) )
		} else if (emis$type == 'norm'){
			isa <- !is.na(x)
			tmp <- apply(matrix(B$L,ncol=J), 2, function(cv) unlist(cov.wt(matrix(x[isa]), cv[isa])[c('cov', 'center')])) # x is already a matrix 
			emis$mu <- tmp['center',]
			emis$var <- tmp['cov', ]				
		}
	}
	return(emis)
}	

splitFarNeighbouryhat<-function(yhat, xPos=NULL, xRange=NULL, maxgap=NULL, state=NULL){
	require(IRanges)
	yhatrle<-Rle(yhat)
	rv<-runValue(yhatrle)
	rl<-runLength(yhatrle)
	ri<-which(rv==as.character(state))
	intStart<-sapply(ri, function(z) sum(rl[seq_len(z-1)])+1)
	intEnd<-sapply(ri, function(z) sum(rl[seq_len(z)]))

	if(length(intStart)>0 && !is.null(maxgap) && !is.null(xPos) || !is.null(xRange)){
		tmp<-splitFarNeighbour(intStart=intStart, intEnd=intEnd, xpos=xPos, xrange=xRange, maxgap=maxgap)
		intStart<-tmp$IS
		intEnd<-tmp$IE		
	}
	return(list(IS=intStart, IE=intEnd))
}


