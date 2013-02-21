\name{biomvRhsmmRun}
\alias{biomvRhsmmRun}
\title{
 Estimating the most likely state sequence using Hidden Semi Markov Model
}
\description{
 This is the working horse of the HSMM run
}
\usage{
biomvRhsmmRun(x, xid = "sampleid", xRange, soj, emis.type = "norm", q.alpha = 0.05, r.var = 0.75, maxit = 1, maxgap = Inf, tol = 1e-06, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	input data matrix or vector
}
  \item{xid}{
	name of the sample
}
  \item{xRange}{
	a \code{IRanges}/\code{GRanges} obejct, same length as \code{x} rows
}
  \item{soj}{
	a list object containing the relevent sojourn distribution parameters
}
  \item{emis.type}{
	type of the emission distribution, only the following types are supported: 'norm', 'mvnorm', 'pois'
}
  \item{q.alpha}{
	a quantile factor controlling the estimated prior for the mean of the emission of each states, \code{seq(from=q.alpha, to=1-q.alpha, length.out=J)}
}
  \item{r.var}{
	a ratio factor controlling the estimated prior for the variance / covariance structure of each states. 
	A value larger than 1 tend to allow larger variation in extreme states;  a value smaller than 1 will decrease the probablity of having extreme state
}
 \item{maxit}{
	max iteration of the EM run with Forward-Backward algorithm
}
  \item{maxgap}{
	max distance between neighbouring feature to consider a split
}
  \item{tol}{
	tolerance level of the likelihood change to determining the termination of the EM run
}
  \item{na.rm}{
	\code{TRUE} if \code{NA} value should be ignored 
}
}
\details{
	from ApplicationNote and Poster
}
\value{
	a list object, 
    \item{\code{yhat}}{the most likely state sequence, same length as \code{x} rows}
    \item{\code{res}}{Object of class \code{"GRanges"} , each range represent one continous segment identified, with sample name slot 'SAMPLE', estimated state slot 'STATE' and segment mean slot 'MEAN' stored in the meta columns}
}
\references{
	Guedon, Y. (2003). Estimating hidden semi-Markov chains from discrete sequences. Journal of Computational and Graphical Statistics, 12(3), 604-639.
}

\seealso{
\code{\link{biomvRseg}}
\code{\link{biomvRseg}}
}

\examples{
	data(coriell)
	x<-coriell[coriell[,2]==1,]
	xgr<-GRanges(seqnames=paste('chr', x[,2], sep=''), IRanges(start=x[,3], width=1, names=x[,1]))
	values(xgr)<-DataFrame(x[,4:5], row.names=NULL)
	xgr<-xgr[order(xgr)]

	J<-2 ; maxk<-50
	soj<-list(J=J, maxk=maxk, type='nbinom', d=cbind(dunif(1:maxk, 1, maxk), dunif(1:maxk, 1, maxk)))
	soj$D <- sapply(1:J, function(j) rev(cumsum(rev(soj$d[1:maxk,j]))))
	sample<-colnames(coriell)[4]
	runout<-biomvRhsmmRun(matrix(values(xgr)[,sample]), sample, xgr, soj, emis.type='norm')
}

\keyword{ hsmm }
