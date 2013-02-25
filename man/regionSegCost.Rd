\name{regionSegCost}
\alias{regionSegCost}
\title{
	 Regional segmentation cost 
}
\description{
	To calculate regional cost matrix for the initial stage and second merging statge of the segmentation model.
}
\usage{
regionSegCost(x, maxk = NULL, segs = NULL, family = NULL, alpha = NULL, useSum = TRUE, useMC = FALSE, comVar = TRUE)
}

\arguments{
  \item{x}{
	The input data matrix or vector
}
  \item{maxk}{
	Maximum number of index to search forward
}
  \item{segs}{
	Starting indices (excluding 1) for the candidiate segments, for the second stage model, \code{maxk} will be overridden with \code{length(segs)+1}.
}
  \item{family}{
	which expontial family the data belongs to, possible values are  'norm', 'pois' and 'nbinom'
}
  \item{alpha}{
	alpha matrix for negative bionomial cost calculation, estimated from \code{\link{regionSegAlphaNB}}
}
  \item{useSum}{
	TRUE if using grand sum across sample / x columns, like in the \code{tilingArray} solution
}
  \item{useMC}{
	TRUE if \code{\link{mclapply}} should be used to speed up
}
  \item{comVar}{
	TRUE if assuming common variance across samples (x columns)
}
}
\details{
	Preparing the cost matrix for the follow-up segmentation. Using residual sum of squares for 'norm' data, and negative log-likehood for 'pois' and 'nbinom' data.
}
\value{
	Matrix with \code{maxk} rows and \code{nrow(x)} columns, or a \code{length(segs)+1} square matrix for the second stage model.
}


\seealso{
	\code{\link{regionSegAlphaNB}}
}
\examples{
	x<-matrix(rnorm(120), ncol=3)
	Ca<-regionSegCost(x, maxk=20, family='norm')
	dim(Ca) # [1] 20 40
	Cb<-regionSegCost(x, segs=as.integer(c(3, 6, 12, 30)), family='norm')
	dim(Cb) # [1] 5 5
}
