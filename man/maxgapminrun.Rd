\name{maxgapminrun}
\alias{maxgapminrun}
\title{
	Max-gap-min-run algorithm for 2 states segmentation
}
\description{
	A custom Max-gap-min-run implementation using physical postion for gap and run length calculation.
}
\usage{
maxgapminrun(x, xpos = NULL, xrange = NULL, cutoff = NULL, q = 0.9, minrun = 5, maxgap = 2, splitLen = Inf)
}

\arguments{
  \item{x}{
	a numeric vector for the input signal
}
  \item{xpos}{
	a numeric vector, same length as x, carrying postional information for each element of x 
}
  \item{xrange}{
	an \code{IRanges} object,  same length as x, carrying range information for each element of x 
}
  \item{cutoff}{
	numeric value used as cutoff, optional if \code{q} is specified
}
  \item{q}{
	numeric value used to derive cutoff of x, as the \code{q} quantile of x , optional if \code{cutoff} is specified
}
  \item{minrun}{
  	minimum run length for the resulting segments
}
  \item{maxgap}{
	numeric value, maximum gap length of region with signal below the cutoff.
}
  \item{splitLen}{
	numeric value, maximum length of segments, split if too long
}
}
\details{
	A custom implementation of max-gap-min-run algorithm.
}
\value{
	a list of segment starts and ends indices
  \item{IS}{the start index for each segment}
  \item{IE}{the end index for each segment}
  \item{CUTOFF}{the cutoff value used in the run}
  \item{MG}{the parameter value for \code{maxgap}}  
  \item{MR}{the parameter value for \code{minrun}}
  \item{SPL}{the parameter value for \code{splitLen}}
}

\author{
Yang Du
}

\seealso{
\code{\link{biomvRhsmm}}
\code{\link{biomvRseg}}
}

\examples{
	x<-rpois(50, 10)
	xpos<-rnorm(50, 300, 100)
	xpos<-xpos[order(xpos)]
	maxgapminrun(x, xpos, cutoff=9.5, maxgap=30, minrun=100)
}

