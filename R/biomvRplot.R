
biomvRplot<-function(x, xMeta=NULL, xAnno=NULL, xGrp=NULL, xSeg=NULL, ftcol=c('green', 'cyan', 'tomato','purple','gold', 'violet'), main='bioMvCNS plot', ylab=expression(paste('Normalized ',italic("log"[2]*""), ' intensity', sep='')), width=12, height=10, ...){
	# core function to plot a region, will provide a higher level wrapper to handle different chr and strand
	# showing nc rows of profile, 1 row of annotation, and 1 for legend
	# x, a matrix or a multivaRseg object,  same defination as for segmentation method
	# xMeta, a GRange object, containing coordinates for each probe, and optional annotation in the elementMetadata
	# 		containg at least 1 column and be the first column for the type of each feature row, 'CDS', 'UTR', 'Intron', 'intergenic'

	# in case there are more than one group, a better strategy would be to use one column for each group xMeta for each column and xAnno for the esimated segments
	# add more rows for each groups in the annotation region, as the predicted outcome
	
	
	# xAnno, a GRange object, containing the annotation for this region, containg at least 1 column and be the first column 
	# 		in the elementMetadata for the feature type, 'CDS', 'UTR', 'Intron', 'intergenic'
	
	## define constants
	if(!is.matrix(x)){
		warning('x is not a matrix, coercing to matrix with 1 column !!!')
	 	x <- matrix(x, ncol=1)
	} 
	nc<-ncol(x)		# number of series 
	nr<-nrow(x)	#number of features

	# check groups
	if(is.null(xGrp)){
		xGrp<-rep(1, nc)
	} else if(is.vector(xGrp) && length(xGrp) != nc){
		stop("length of vecter xGrp should be the same as ncol(x) !!!")
	}
	ng<-length(unique(xGrp))
	xfts<-fts<-character(0)
	xhascol<-fthascol<-FALSE
	## somechecking or range object
	
	if( !is.null(xMeta) && length(xMeta)==nr && ncol(values(xMeta)) >= ng && 
	all(sapply(1:ng, function(g) is.character(values(xMeta)[,g]))) && 
	length(unique(unlist(lapply(1:ng, function(g) unique(values(xMeta)[,g])))))<=5){
		xhascol<-TRUE
		xfts<-unique(unlist(lapply(1:ng, function(g) unique(values(xMeta)[,g]))))
		if(!is.null(xAnno)){
			if(class(xAnno)=='GRanges' && ncol(values(xAnno)) >0  && 
			is.character(values(xAnno)[,1]) && length(unique(values(xAnno)[,1]))<=5 ){
				xAnno<-xAnno[start(xAnno)>=min(start(xMeta)) & end(xAnno)<=max(end(xMeta))]
				fthascol<-TRUE
				fts<-c(fts, unique(values(xAnno)[,1]))
			} else if (class(xAnno)=='GRangesList' && length(xAnno) == ng &&
			all(sapply(xAnno, function(z) ncol(values(z))>0 && is.character(values(z)[,1]) && 
			length(unique(unlist(lapply(xAnno, function(z) values(z)[,1]))))<=5) ) ){
				fthascol<-TRUE
				fts<-c(fts, unique(unlist(lapply(xAnno, function(z) values(z)[,1]))))
				xAnno<-lapply(xInfo$xAnno, function(z) z[start(z)>=min(start(xInfo$xMeta)) & end(z)<=max(end(xInfo$xMeta))])
			} else {
				warning('the class of xAnno is not valid, ignored !!')
			}
		}
	} else {
		warning('no valid xMeta found, plot features sequentially, xAnno will be ignored if specified !!!')
		xMeta<-GRanges(seqnames = Rle('dummy', nr), ranges = IRanges(start=seq_len(nr), width=1), strand =Rle('*', nr) )
		if(!is.null(xAnno)){
			xAnno<-NULL
		}
	}

	
	## plot limits x
#	xlim<-c(floor(min(start(xMeta), if(is.null(xAnno)) NA else start(xAnno), na.rm=T)/10)*10, ceiling(max(end(xMeta), if(is.null(xAnno)) NA else end(xAnno), na.rm=T)/10)*10)
	xlim<-c(floor(min(start(xMeta), na.rm=T)/10)*10, ceiling(max(end(xMeta), na.rm=T)/10)*10)
	## feature color coding
	if(xhascol || fthascol){
		typecode<-unique(c(xfts, fts))
		typecode<-typecode[order(typecode)]
		colcode<-ftcol[seq(along=typecode)]
		names(colcode)<- typecode
	} else {
		 colcode<-ftcol[0]
		 names(colcode)<-character(0)
	}
	
	# plot device setup
	if(is.null(dev.list())){
		dev.new(width=width, height=height)
	} else {
		if(dev.size()[2]<12) warning('You may need a larger/taller device to plot on.')
	}
	har<-c(rep(5/6/nc, nc), rep(1/12, 2)) # the aspect ratio of regions heights
	war<-c(1/12,11/12)	# the aspect ratio of regions width
	nf<-layout(matrix(c(rep(nc+3,nc), rep(nc+4, 2), seq_len(nc+2)), nc+2, 2), heights=har, width=war) 
	par(mar=c(3,2.5,1,2.5),oma=c(0.5,0.5,4.5,0.5)) # the margin may need to be changed
	#par(mar=c(0.5,0.5,0.5,0.5),oma=c(0.5,0.5,4.5,0.5)) # the margin may need to be changed
	#layout.show(nf) 
	
	# plot profile @ 1:nc
	for (i in 1:nc){
		ylim<-c(floor(min(x[,i]))-1, ceiling(max(x[,i]))+1) 
		## use the same color coding for annotation if avaliable
		if(xhascol){
			# plot the first one like NA ones
			if(!is.na(xfts[1])){
				coli<-values(xMeta)[,xGrp[i]]==xfts[1]
				plot.default(y=x[coli,i], x=(start(xMeta[coli,])+end(xMeta[coli,]))/2, xlim=xlim, ylim=ylim, axes=FALSE, main=NULL, xlab='',ylab='', col=as.character(colcode[xfts[1]]))		
			} else{
				coli<-is.na(values(xMeta)[,1])
				plot.default(y=x[coli,i], x=(start(xMeta[coli,])+end(xMeta[coli,]))/2, xlim=xlim, ylim=ylim, axes=FALSE, main=NULL, xlab='',ylab='', col=as.character(colcode[is.na(names(colcode))]))		
			}

			if(length(xfts)>1){
				for(fti in 2:length(xfts)){
					coli<-which(values(xMeta)[,xGrp[i]]==xfts[fti])
					points(y=x[coli,i], x=(start(xMeta[coli,])+end(xMeta[coli,]))/2, col=as.character(colcode[xfts[fti]]))
				}
			}
			#point the rest
		} else {
			plot.default(y=x[,i], x=(start(xMeta)+end(xMeta))/2, xlim=xlim, ylim=ylim, axes=FALSE, main=NULL, xlab='',ylab='')
		}
		axis(side=2, at=seq(ylim[1], ylim[2], length.out=5), labels=seq(ylim[1], ylim[2], length.out=5), tick=TRUE, cex.axis = 2)	
		axis(side=4, at=mean(ylim[1], ylim[2]), labels=paste('s',i,' / g',xGrp[i],sep=''), tick=FALSE, cex.axis = 2, padj=0)	
		# segmentation index of x transform to xMeta
		for(xi in xSeg[[i]]){
			lines(rep(start(xMeta[xi]),2), ylim, col='grey')
		}
		#abline(v=start(xMeta[xSeg[[i]]]), col='grey')
	}
	
	#plot  annotation @ nc+1 if there is any
	# stack  anno by group
	aH<-2
	plot.default(x=xlim,y=c(0,0),xlim=xlim, ylim=c(0,aH),axes=FALSE,main=NULL,type='l', xlab='',ylab='')
	axis(side=1, at=seq(xlim[1],xlim[2], length.out=10), labels=as.integer(seq(xlim[1],xlim[2], length.out=10)), tick=TRUE, cex.axis = 2)	
	if(fthascol){
		if(class(xAnno)=='GRanges'){
			rect(start(xAnno), 0, end(xAnno), aH, col=as.character(colcode[values(xAnno)[, 1]]), border = NA)
		} else{
			# stack up
			for(g in 1:ng){
				rect(start(xAnno[[g]]), aH-(g-1)*aH/ng, end(xAnno[[g]]), aH-g*aH/ng, col=as.character(colcode[values(xAnno[[g]])[, 1]]), border = NA)
			}	
		}
		axis(side=2, at=(seq_len(ng)-0.5)*aH/ng, labels=paste('g', seq_len(ng),sep=''), tick=F, cex.axis = 1.5, las=1)	
	}

	
	#plot legend @ nc+2
	plot.default(xlim,c(0,0),xlim=xlim, ylim=c(-1,1),axes=FALSE,main=NULL,type='n', xlab='',ylab='')
	nleg<-length(colcode)
	if(nleg>0){
		legpos<-seq(xlim[1], xlim[2], length.out=nleg+1)
		legcentpos<-(legpos[1:(nleg)]+legpos[2:(nleg+1)])/2
		legdist<-mean(legpos[2:(nleg+1)]-legpos[1:(nleg)])
		rect(legcentpos-legdist/6, rep(-1,nleg),legcentpos, rep(1,nleg), col=as.character(colcode), border = NA)
		text(x=legcentpos+legdist/6,y=0,label=names(colcode), cex=2.5)
	}


	#plot profile label @ nc+3
	plot.default(0,0,axes=F,xlab='', ylab='', type='n')
	text(0,0, ylab, srt=90, cex=3) 
	
	#plot legend label @ nc+4
	plot.default(0,0,axes=F,xlab='', ylab='', type='n')
	text(0,0,'Pos.', cex=1.5)
	## add a main title at top
	mtext(main, side=3, line=1, cex=3, outer=TRUE) 
}




biomvRGviz<-function(exprgr, gmgr=NULL, prange=NULL, regionID=NULL, seggr=NULL, plotstrand='*', eps=TRUE, showId=TRUE, ...){
	# exprgr, probe info, with first data column as expression value
	# gmgr,  related annotation data, optional, a TYPE mcol has to be there...fragile
	# seggr, segmentation info, optional, a STATE mcol has to be there...fragile
	# prange, a range to plot, optional, a reasonable region is strongly advisable.
	#Ideogram track obviously is not avaliable for pig on UCSC
	#ideoTrack <- IdeogramTrack(genome = "susScr3", chromosome = "chrX")
	options(ucscChromosomeNames=FALSE)
	#check if there are more chrs, and no prange
	if(length(unique(as.character(seqnames(exprgr))))>1 && is.null(prange)){
		stop("More than 1 chr in the Grange object, yet no plot region defined! ")
	}
	if(is.null(prange)){
		prange<-c(unique(as.character(seqnames(exprgr))), floor(min(start(exprgr))/1000)*1000, ceiling(max(end(exprgr))/1000)*1000)
	}
	# handle the color automatically according to gmgr and seggr
	colors<-c('cyan', 'tomato', 'green','purple','gold', 'violet')
	typecode<-NULL
	if(!is.null(gmgr)){
		typecode<-c(typecode,unique(values(gmgr)[,'TYPE']))
	}
	if(!is.null(seggr)){
		typecode<-unique(c(typecode, unique(values(seggr)[,'STATE'])))
	}
	if(!is.null(typecode)){
		typecode<-typecode[order(typecode)]
		params <- as.list(colors[seq_along(typecode)])
		names(params)<-typecode
	} else {
		params<-list()
	}
	
	if(hasArg(ylab))  ylab <- list(...)$ylab else ylab<-NULL
	if(hasArg(main))  main <- list(...)$main else main<-NULL
	
	if(is.null(main))	main<-paste(ifelse(is.null(regionID), 'region', gene), '@', prange[1], '.', prange[2],'-', prange[3], '@', colnames(mcols(exprgr)),  sep='')
	trackList<-list()
		
	# datatrack + 
	if(plotstrand == '+' | plotstrand == '*'){
		if(is.null(ylab)) ylab<-paste('E@+', sep='') else ylab<-paste(ylab, '@+', sep='')
		dpTrack <- DataTrack(exprgr[seqnames(exprgr)==prange[1] & (strand(exprgr)=='+' | strand(exprgr)=='*') & start(exprgr) >= as.numeric(prange[2]) & end(exprgr) <= as.numeric(prange[3])],  name = ylab, background.title = "darkblue", type = c('p'),  legend = TRUE)				
		trackList<-append(trackList, dpTrack)
		if(!is.null(seggr)){
			# segmentation + as a seprate state annotation
			# no id, with legend
			segp<-seggr[seqnames(seggr)==prange[1] & (strand(seggr)=='+' | strand(seggr)=='*') & start(seggr) >= as.numeric(prange[2]) & end(seggr) <= as.numeric(prange[3])]
			spTrack<- AnnotationTrack(segp,  name=paste('S', '@+', sep=''), id=values(segp)[,'STATE'] ,background.title = "Gray", background.panel = "#FFFFFF", showFeatureId = showId, shape = "box")
			feature(spTrack)<- values(segp)[,'STATE']
			trackList<-append(trackList, spTrack)
		}
		
		if(!is.null(gmgr)){
			# annodat track +
			gmp<-gmgr[seqnames(gmgr)==prange[1] & (strand(gmgr)=='+' | strand(gmgr)=='*') & start(gmgr) >= as.numeric(prange[2]) & end(gmgr) <= as.numeric(prange[3])]
			apTrack<- AnnotationTrack(gmp, group=names(gmp), name=ifelse(length(gmp)==0, '', paste('G', '@+', sep='')), background.title = "brown", background.panel = "#FFFEDB", showId = showId, shape = "box")
			feature(apTrack)<- values(gmp)[,'TYPE'] ## this is a strong requirement
			trackList<-append(trackList, apTrack)
			rm(gmp)
		}				
	}
	
	# axis track
	axisTrack <- GenomeAxisTrack(add53 = TRUE, add35 = TRUE, littleTicks = TRUE)
	trackList<-append(trackList, axisTrack)
	
	if(plotstrand == '-' | plotstrand == '*'){			
		if(!is.null(gmgr)){
			# annodat track -
			gmm<-gmgr[seqnames(gmgr)==prange[1] & (strand(gmgr)=='-' | strand(gmgr)=='*') & start(gmgr) >= as.numeric(prange[2]) & end(gmgr) <= as.numeric(prange[3])]
			amTrack <- AnnotationTrack(gmm, group=names(gmm), name=ifelse(length(gmm)==0, '', paste('G', '@-', sep='')), background.title = "brown", background.panel = "#FFFEDB", showId = showId, shape = "box")
			feature(amTrack)<- values(gmm)[,'TYPE']
			trackList<-append(trackList, amTrack)
			rm(gmm)
		}
		if(!is.null(seggr)){
			# segmentation + as a seprate state annotation
			# no id, with legend
			segm<-seggr[seqnames(seggr)==prange[1] &  (strand(seggr)=='-' | strand(seggr)=='*')  & start(seggr) >= as.numeric(prange[2]) & end(seggr) <= as.numeric(prange[3])]
			smTrack<- AnnotationTrack(segm,  name=paste('S', '@-', sep=''),  id=values(segm)[,'STATE'], background.title = "Gray", background.panel = "#FFFFFF", showFeatureId = showId, shape = "box")
			feature(smTrack)<- values(segm)[,'STATE']
			trackList<-append(trackList, smTrack)
		}
		#datatrack - 
		if(is.null(ylab)) ylab<-paste('E@-', sep='') else ylab<-paste(ylab, '@-', sep='')
		dmTrack <- DataTrack(exprgr[seqnames(exprgr)==prange[1] & (strand(exprgr)=='-' | strand(exprgr)=='*') & start(exprgr) >= as.numeric(prange[2]) & end(exprgr) <= as.numeric(prange[3])], name = ylab, background.title = "darkblue", type = c('p'),   legend = TRUE)
		trackList<-append(trackList, dmTrack)	
	}

	# start ploting
	# initial grapic dev.
	graphics.off()
	if(eps){
		setEPS()
		postscript(paste(main, '.', plotstrand,'.eps', sep=''), paper='special', width=25, height=10, horizontal=F, fonts=c("sans"))
	} else {
		pdf(paste(main, '.', plotstrand, '.pdf', sep=''), width=25, height=10)
	}
	# append general plotting params.
	params <- append(params, list(fontsize=20, cex.title=1.5, main=main))
	do.call(plotTracks,c(list(trackList), params)) 
	dev.off()
}

