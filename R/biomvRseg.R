## main function for normized numeric intensity vectors
biomvRseg<-function(x, maxk=NULL, maxseg=NULL, family='norm', grp=NULL, penalty='BIC', clusterm=NULL, twoStep=TRUE, segDisp=FALSE, useMC=FALSE, useSum=TRUE, comVar=TRUE){
	# input to the main function is a matrix object, x, features on the same strand and same chr
	# optional grouping factor, length of which should be the same as the column of x
	# input checking and preparion
	if (is.null(family) || (family != 'norm' && family != 'pois' && family != 'nbinom')) 
    	stop("'family' must be specified, currently only 'norm' ,'pois' and 'nbinom' are supported !")
	
	if (!is.numeric(x) || !(is.vector(x) || is.matrix(x)) || class(x)=='GRanges') 
        stop("'x' must be a numeric vector or matrix or a GRanges object.")
    if(class(x)=='GRanges') {
    	xid<-names(values(x))
    	xRange<-ranges(x)
    	x<-as.matrix(values(x))
    }    
	if(length(dim(x))==2){
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
	
	if (is.null(maxseg) || !is.numeric(maxseg) || (length(maxseg) != 1) || (maxk <= 1) ||  (maxk > nr)) 
		 stop(sprintf("'maxseg' must be a single integer between 2 and the number of rows of 'x': %d.", nr))
	maxseg<-as.integer(maxseg) # this is the number for initial candidate region
 	
	 if (is.null(maxk) || !is.numeric(maxk) || (length(maxk) != 1) || (maxk <= 1) ||  (maxk > nr)) 
	 	 stop(sprintf("'maxk' must be a single integer between 2 and the number of rows of 'x': %d.", nr))
	maxk<-as.integer(maxk) # this is the maximal windows size for the initial segmentation
	
	if(maxk*maxseg<nr) 
		stop(sprintf("the product of 'maxseg' and 'maxk' is smaller than the number of rows of 'x': %d.", nr))

    # penalty
    penaltymethods<-c('none','AIC','AICc','BIC','SIC','HQIC')
    penalty<-match.arg(penalty, penaltymethods)
    if(penalty=='SIC') penalty='BIC' 
    
    # check grp setting, cluster if needed, otherwise treat as one group	
   	if(!is.null(grp)) grp<-as.character(grp)
	grp<-preClustGrp(x, grp=grp, clusterm=clusterm)
	
	## initialize the output vectors
	segStart <- vector(mode="list", length=nc)
	segMean <- vector(mode="list", length=nc)
	## start the initial segmetation within each grp
	for(g in unique(grp)){
		cat(sprintf("step 1 segmentation for group %s\n", g))
		d<-sum(grp==g)
		## construct costmatrix, family wise
		if(family=='nbinom'){
			if (segDisp){
				alpha<-regionSegAlphaNB(x[,which(grp==g)], maxk=maxk, useMC=useMC)
			} else {
				dispersion<-estimateSegCommonDisp(x[,which(grp==g)])
				alpha<-matrix(dispersion, maxk, nr)
			}
			C<-regionSegCost(x[,which(grp==g)], maxk=maxk, family=family, alpha=alpha, useSum=useSum, useMC=useMC)
		} else {
			C<-regionSegCost(x[,which(grp==g)], maxk=maxk, family=family, useSum=useSum, comVar=comVar)
		}
#		Res<-.Call("univaRseg", C, maxseg,  PACKAGE = "biomvCNS")
		Res<- .C("univaRseg", as.double(C), as.integer(maxseg), as.integer(maxk), as.integer(nr), cost=double(maxseg*nr), pos=integer((maxseg-1)*nr), minC=double(maxseg), segS=integer(maxseg*maxseg))#, PACKAGE = "biomvRCNS")
		#Res$minC, mincost for each number of changepoint option
		#Res$segS, the start position of each segment, 1 excluded and n+1 asserted

		if( family== 'norm'){
			# norm
			#logL<- -nr/2*(1+log(2*pi)+log(Res$minC/nr)) # in tilingArray is a mixture of common and change variance
			if(comVar){
				logL<- -nr/2*(d+d*log(2*pi)+log(Res$minC/nr)) # the likelihood ratio cost, p104
			} else {
				logL<- -(nr+nr*d*log(2*pi)+Res$minC)/2 # the likelihood ratio cost, p134
			}
		} else {
			# pois and nb
			logL<- 0-Res$minC
		}
		
		# number of parameters are different, pois has only 1 for each segment, whereas norm and NB both have two
		# for norm, comVar=T, useSum=T, (s-1)*pos + s*mean + 1*sd = 2s
		# for norm, comVar=F, useSum=T, (s-1)*pos + s*mean + s*sd = 3s-1
		# for norm, comVar=T, useSum=F, (s-1)*pos + s*d*mean + 1*d*sd = (s+1)(d+1)-2
		# for norm, comVar=F, useSum=F, (s-1)*pos + s*d*mean + s*d*sd = 2*s*d+s-1 
		# for pois, (s-1)*pos + s*mean=2s-1
		# for NB, segDisp=T, (s-1)*pos + s*mean + s*alpha = 3s-1
		# for NB, segDisp=F, (s-1)*pos + s*mean + 1*alpha = 2s
		
		if(family=='pois'){
			nP<-seq_len(maxseg)*2-1
		} else if ( (family=='nbinom' && !segDisp) || (family=='norm' && comVar && useSum)){
			nP<-seq_len(maxseg)*2
		} else if ( (family=='nbinom' && segDisp) || (family=='norm' && !comVar && useSum)) {
			nP<-seq_len(maxseg)*3-1
		} else if ( family=='norm' && comVar && !useSum) {
			nP<-(seq_len(maxseg)+1)*(d+1)-2
		} else if ( family=='norm' && !comVar && !useSum) {
			nP<-seq_len(maxseg)*d*2+seq_len(maxseg)-1
		}
		
#	BICpenalty = -0.5*sum(log(nTotalWindow)) - nCpts*log(length(combX)) + 0.5*log(nTotal)
#	mBIC = lik1-lik0+BICpenalty

		
		# the optimal number of final segments for series within this group
		rN<-switch(penalty,
			none = which.max(logL),
			AIC = which.min(sapply(seq_len(maxseg), function(x) -2*logL[x]+2*nP[x])),
			AICc = which.min(sapply(seq_len(maxseg), function(x) -2*logL[x]+2*nP[x]*(nP[x]+1)/(nr*d-nP[x]-1))),
			BIC = which.min(sapply(seq_len(maxseg), function(x) -2*logL[x]+nP[x]*log(nr*d))),
			HQIC = which.min(sapply(seq_len(maxseg), function(x) -2*logL[x]+2*nP[x]*log(log(nr*d)))),
			stop('Invalid value argument for penalty'))	
		cat(sprintf("Step 1 segmentation complete for group %s\n", g))
		
		#rN could be supplied by the user, thus only keeping the initial segment candidate
		
		# now apply individually for each serie within this group
		for(c in which(grp==g)){
			## here add checking whether a 2nd step is needed
			if(d==1 || !twoStep || rN==maxseg){
				cat(sprintf("No need to run step 2 merging for column %s from group %s\n", c, g))
				## no need to do 2nd step
				rIdx<-c(mat2list(Res$segS, maxseg)[[rN]])
				j<-c(1,rIdx) # col, segStart
				i<-c(rIdx-1, nr) # segEnd
				k<-  i - j+1 # row for C and alpha matrix
				
				segStart[[c]]<-rIdx
				segMean[[c]]<-sapply(seq_len(rN), function(z) mean(x[i[z]:j[z],c]))
				cat(sprintf("Processing complete for column %s from group %s\n", c, g))
			} else {
				cat(sprintf("Step 2 merging for column %s from group %s\n", c, g))
				# to run a 2nd step merging
				# get all candidates coordinates
				rRegs<-mat2list(Res$segS, maxseg)[[maxseg]]	
				# generate a regional cost for each series
				if(family=='nbinom'){
					if(segDisp){
						ralpha<-regionSegAlphaNB(x[,c], segs=rRegs, useMC=useMC)
					} else {
						ralpha=matrix(dispersion,maxseg, maxseg)
					}
					rC<-regionSegCost(x[,c], segs=rRegs, family=family, alpha=ralpha, useSum=useSum, useMC=useMC) 
				} else  {
					rC<-regionSegCost(x[,c], segs=rRegs, family=family, useSum=T, comVar=comVar)
				}		
				# use this regional cost to do dp merging
#				rRes<-.Call("univaRseg", rC , as.integer(rN),  PACKAGE = "biomvRCNS")
				rRes<- .C("univaRseg", as.double(rC), as.integer(rN), as.integer(maxseg), as.integer(maxseg), cost=double(rN*maxseg), pos=integer((rN-1)*maxseg),minC=double(rN*maxseg), segS=integer(rN*rN))#, PACKAGE = "biomvRCNS")
				
				rIdx<-mat2list(rRes$segS, rN)[[rN]]
				j<-c(1,rRegs[rIdx-1]) # col, segStart, original index, including 1
				i<-c(rRegs[rIdx-1]-1, nr) # segEnd, original index
				
				# the result here need to be put into the final returning value
				segStart[[c]]<-rRegs[rIdx-1]				
				segMean[[c]]<-sapply(seq_len(rN), function(z) mean(x[i[z]:j[z],c]))	
				cat(sprintf("Step 2 merging complete for column %s from group %s\n", c, g))
			} # end 2nd step if
		} # end c for
	} # end g for
	new("biomvRseg",
    x = x,
    segStart = segStart,
    segMean = segMean,
    group=as.integer(grp),
	family=family)
}



##################################################
# segmentation utility functions
##################################################

regionSegCost<-function(x, maxk=NULL, segs=NULL, family=NULL, alpha=NULL, useSum=TRUE, useMC=FALSE, comVar=TRUE, naVal=.Machine$double.xmax){
	# the same cost function for both initial DP segmentation and 2nd DP merging
	# segs, a vector of the starting index of each candiate region, except for the first position 1
	# check input x, convert to vecter if necessary
	if (is.null(family) || (family != 'norm' && family != 'pois' && family != 'nbinom')) 
        stop("'family' must be specified, currently only 'norm' ,'pois' and 'nbinom' are supported !")
  
	if (!is.numeric(x) || !(is.vector(x) || is.matrix(x))) 
        stop("'x' must be a numeric vector or matrix.")
    if (is.vector(x)) {
        d <- 1
        r <- x
        if( family== 'norm'){
        	 q<- x*x
        }
        x<-matrix(x, ncol=d)
    } else {
    	d <- ncol(x)
        r <- rowSums(x)
        if( family== 'norm'){
       	 q <-rowSums(x*x)
        }
    }
    n<-length(r)
    
    ## add a checking for useSum in norm models
    if(d==1 && !useSum){
    	useSum<-TRUE
    	warning("For univariate normal data, 'useSum' reset to TRUE, though result should be identical.")
    }
    
    ## check input parameter maxk and segs
     if (is.null(segs)) {
     	if (is.null(maxk)){
     		# both not specified, offer a warning
     		warning(sprintf("both 'maxk' and 'segs' are not specified, a full size %d x %d cost matrix would be generated.", n, n))
			maxk<-n
     	} else {
     		# cost for the fisrt DP seg step
			 if (!is.numeric(maxk) || (length(maxk) != 1) || (maxk <= 1) ||  (maxk > n)) 
			    stop(sprintf("'maxk' must be a single integer between 2 and the number of rows of 'x': %d.", n))
     	}
		segs<-seq_len(n-1)+1
		N<-n
     } else if ( is.numeric(segs) && min(segs)>1 && max(segs)<=n  ) {
     	## if no segs specified, treat every position individually, no grouping, defult for the 1st step
    	N<-length(segs)+1
     	if (!is.null(maxk)){
     		# both specified, offer a warning
     		warning(sprintf("both 'maxk' and 'segs' are specified, maxk would be overrided with lenth(segs)+1=%d", N))
     	}
     	maxk<-N 
     } else {
     	stop(sprintf("'segs' must be a vector of integers between 2 and the number of rows of 'x': %d.", n))
     }
     
    ## check nb and alpha,  alpha has to be consistent with C
    if( family == 'nbinom' && (!is.matrix(alpha) || nrow(alpha) != maxk && ncol(alpha)!=N))
		stop(sprintf("'alpha' must be a matrix of %d x %d .", maxk, N))
     
    # all end position of each candidate region
	e<-c(segs-1, n)
    
    # sum(x) wrt segs, for tilingArray and likelihood ratio cost 1
	cr<-cumsum(r)
	crs<-rep(cr[e[1]],N)
	crs[2:N]<-cr[e[2:N]]-cr[e[1:(N-1)]]
	crss<-cumsum(crs)
	
	# generate the cost wrt segs and maxk
	C <- matrix(as.numeric(naVal), nrow = maxk, ncol = N)
	sk<-cumsum(c(segs-1,n)-c(1,segs)+1)
	k <- 1:maxk 
	
    if( family== 'norm'){
    	# sum(x^2) wrt segs for norm only
		cq<-cumsum(q)
		cqs<-rep(cq[e[1]],N)
		cqs[2:N]<-cq[e[2:N]]-cq[e[1:(N-1)]]
		cqss<-cumsum(cqs)
		
		## there is a difference between grand mean and mean vector, currently, we are using grand mean
		## a complete solution would be a sum vector within a sum matrix, and then sum its segments up
		if(useSum){
			if(comVar){
				C[, 1] <- cqss[k] - crss[k] * crss[k]/(sk[k]*d) # old tilingArray solution, assuming common variance and use grand sum
			} else {
				C[, 1] <- log((cqss[k] - crss[k] * crss[k]/(sk[k]*d))/sk[k])*sk[k] # likelihood ratio cost, assuming different variance and use grand sum
			}
		} else {
			## sum(x) wrt segs, for the new likelihood ration cost
			cx<-sapply(1:d, function(j) cumsum(x[,j]))
			cxs<-matrix(, N, d)
			cxs[1, ]<- cx[e[1], ]
			cxs[2:N, ]<- sapply(1:d, function(j) cx[e[2:N],j]- cx[e[1:(N-1)],j])
			cxss<-sapply(1:d, function(j) cumsum(cxs[,j]))
		
			if(comVar){
				C[, 1] <- cqss[k] - sapply(k, function(i) sum(cxss[i,]^2))/sk[k]	# likelihood ratio cost, column wise cumsum and assume common variance
			} else {
				C[, 1] <- log((cqss[k] - sapply(k, function(i) sum(cxss[i,]^2))/sk[k])/sk[k])*sk[k]	# new likelihood ratio cost, column wise cumsum and different variance
			}
		}
		for (k in 1:(ifelse(N==maxk,maxk-1,maxk)) ) {
		    i <- 1:(N - k)
		    cqssk <- cqss[i + k] - cqss[i]
		    skk<-sk[i+k]-sk[i]
	    	if(useSum){
	    		crssk <- crss[i + k] - crss[i]
	    		if(comVar){
	    			C[k, i+1] <- cqssk - crssk * crssk/(skk*d) # old tilingArray solution, use grand sum, common variace, and grand sum
	    		} else {
	    			C[k, i+1] <-  log((cqssk -  crssk * crssk/(skk*d))/skk)*skk #  likelihood ratio cost, assuming different variance and use grand sum
	    		}
	    	} else {
	    		cxssk<-matrix(t(sapply(i, function(z)  cxss[z + k,] - cxss[z,])), ncol=d, nrow=N-k)
	    		if(comVar){
	    			C[k, i+1] <- cqssk - sapply(i, function(z) sum(cxssk[z,]^2))/skk # likelihood ratio cost, column wise cumsum, common variance
	    		} else {
	    			C[k, i+1] <- log((cqssk - sapply(i, function(z) sum(cxssk[z,]^2))/skk)/skk)*skk # likelihood ratio cost, column wise cumsum, variance change
	    		}
	    	}
		}
		if(N==n) C[1,]<-naVal #  give NA/Inf to the cost of single point in the first step, to avoid a positive epsilon or NaN, thus exclude consecutive changes
		
    } else if (family== 'pois') { # for poisson
    	C[, 1] <- crss[k]*log(crss[k]/(sk[k]*d))
		for (k in 1:(ifelse(N==maxk,maxk-1,maxk)) ) {
		    i <- 1:(N - k)
		    crssk <- crss[i + k] - crss[i]
		    skk<-sk[i+k]-sk[i]
		    C[k, i+1] <- crssk*log(crssk/(skk*d))
		}
		#C[1,]<-NA ## experimental, to avoid consecutive changes, doesnot work well, and inf cause 
		C<- -C
    }  else { # for negative bionomial
    	if(all(alpha==0)){
    		C[, 1] <- sapply(k, function(z) NBlogL1stTerm(c(x[seq_len(e[z]),]), alpha[z,1])) + crss[k]*log(crss[k]/(sk[k]*d)) - log((1+alpha[,1]*crss[k]/(sk[k]*d))^(crss[k]+(sk[k]*d)/alpha[,1]))
    	} else if(useMC) {
    		 C[, 1] <- unlist(mclapply(sapply(k, function(y) list(y)), function(z) NBlogL1stTerm(c(x[seq_len(e[z]),]), alpha[z,1]))) + crss[k]*log(crss[k]/(sk[k]*d)) - (crss[k]+(sk[k]*d)/alpha[,1])*log(1+alpha[,1]*crss[k]/(sk[k]*d))
    	} else {
    		C[, 1] <- sapply(k, function(z) NBlogL1stTerm(c(x[seq_len(e[z]),]), alpha[z,1])) + crss[k]*log(crss[k]/(sk[k]*d)) - (crss[k]+(sk[k]*d)/alpha[,1])*log(1+alpha[,1]*crss[k]/(sk[k]*d))
    	}  	
    	if(!useSum){
    		C[, 1] <- C[, 1]/(sk[k]*d)	
    	}				
		for (k in 1:(ifelse(N==maxk,maxk-1,maxk)) ) {
			i <- 1:(N - k)
			crssk <- crss[i + k] - crss[i]
			skk<-sk[i+k]-sk[i]
			if(all(alpha==0)){
				C[k, i+1] <- sapply(i, function(z) NBlogL1stTerm(c(x[segs[z+k-1]:e[z+k],]), alpha[k, z+1])) + crssk*log(crssk/(skk*d)) - log((1+alpha[k, i+1]*crssk/(skk*d))^(crssk+(skk*d)/alpha[k, i+1]))
			} else if(useMC) {
				C[k, i+1] <- unlist(mclapply(sapply(i, function(y) list(y)), function(z) NBlogL1stTerm(c(x[segs[z+k-1]:e[z+k],]), alpha[k, z+1]))) + crssk*log(crssk/(skk*d)) - (crssk+(skk*d)/alpha[k, i+1])*log(1+alpha[k, i+1]*crssk/(skk*d))
			}	else {
				C[k, i+1] <- sapply(i, function(z) NBlogL1stTerm(c(x[segs[z+k-1]:e[z+k],]), alpha[k, z+1]))+ crssk*log(crssk/(skk*d)) - (crssk+(skk*d)/alpha[k, i+1])*log(1+alpha[k, i+1]*crssk/(skk*d))
			}
			if(!useSum){
				C[k, i+1] <- C[k, i+1]/(skk*d)		
			} 				
		}
		C<- -C
    }
	return(C)
}


NBlogL1stTerm<-function(x, alpha){
	# here alpha is a sinlge value, estimated from the same segment of the x matrix
	# log first and sum latter
	if(alpha==0){
		0
	} else {
		sum(sapply(x, function(y) if(y==0) 1 else sum(sapply(seq_len(y)-1, function(z) log(z*alpha+1)))))
	}
}


regionSegAlphaNB<-function(x, maxk=NULL, segs=NULL, useMC=FALSE){
	# generate alpha matrix for NB model, in both steps of segmentation cost prep.
	# segs, a vector of the starting index of each candiate region, except for the first position 1
	# check input x, convert to vecter if necessary
	if (!is.numeric(x) || !(is.vector(x) || is.matrix(x))) 
        stop("'x' must be a numeric vector or matrix.")
    if (is.vector(x)) {
        d <- 1
        r <- x
        x<-matrix(x, ncol=d)
    } else {
    	d <- ncol(x)
        r <- rowSums(x)
    }
    n<-length(r)
    
    ## check input parameter maxk and segs
     if (is.null(segs)) {
     	if (is.null(maxk)){
     		# both not specified, offer a warning
     		warning("both 'maxk' and 'segs' are not specified, a full size n x n cost matrix would be generated.")
			maxk<-n
     	} else {
     		# cost for the fisrt DP seg step
			 if (!is.numeric(maxk) || (length(maxk) != 1) || (maxk <= 1) ||  (maxk > n)) 
			    stop(sprintf("'maxk' must be a single integer between 2 and the number of rows of 'x': %d.", n))
     	}
		segs<-seq_len(n-1)+1
		N<-n
     } else if ( is.numeric(segs) && min(segs)>1 && max(segs)<=n  ) {
     	## if no segs specified, treat every position individually, no grouping, defult for the 1st step
    	N<-length(segs)+1
     	if (!is.null(maxk)){
     		# both specified, offer a warning
     		warning("both 'maxk' and 'segs' are specified, maxk would be overrided with lenth(segs)+1")
     	}
     	maxk<-N 
     } else {
     	stop(sprintf("'segs' must be a vector of integers between 2 and the number of rows of 'x': %d.", n))
     }
     # all end position of each candidate region
	e<-c(segs-1, n)
	
	# generate alpha wrt segs and maxk
	alpha <- matrix(as.numeric(NA), nrow = maxk, ncol = N)
	k <- 1:maxk 

	if(useMC){
		alpha[, 1] <-unlist(mclapply(sapply(k, function(y) list(y)), function(z) estimateSegCommonDisp(x[seq_len(e[z[1]]),])))	
	} else {
		alpha[, 1] <-sapply(k, function(z) estimateSegCommonDisp(x[seq_len(e[z]),]))		
	}	
	for (k in 1:(ifelse(N==maxk,maxk-1,maxk)) ) {
		## need a provide a vector for alpha, this is for each segment length k,  estimated from the all segments start from i to j	
	    i <- 1:(N - k)  
		if(useMC){
			alpha[k, i+1] <-unlist(mclapply(sapply(i, function(y) list(y)), function(z) estimateSegCommonDisp(x[segs[z[1]+k-1]:e[z[1]+k],])))
		} else {
			alpha[k, i+1]<-sapply(i, function(z) estimateSegCommonDisp(x[segs[z+k-1]:e[z+k], ]))
		}
	}
	return(alpha)
}



estimateSegCommonDisp<-function(xSeg,tol=1e-06){
	# function imported from edgeR
	# xSeg is a matrix object or a vector
		disp <- 0
		for(i in 1:2) {
			delta <- optimize(commonCondLogLikDerDelta, interval=c(1e-4,100/(100+1)), tol=tol, maximum=TRUE, y=list(xSeg), der=0, doSum=FALSE)
			delta <- delta$maximum
			disp <- delta/(1-delta)
		}
		disp
}

