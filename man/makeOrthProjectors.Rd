% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeOrthProjectors.R
\name{makeOrthProjectors}
\alias{makeOrthProjectors}
\title{Construct Orthogonal Projector Matrices}
\usage{
makeOrthProjectors(BlkDesList)
}
\arguments{
\item{BlkDesList}{a list of block design matrices generated by
\code{\link{makeBlkDesMat}}.}
}
\value{
A list containing matrices.
}
\description{
Construct a list of orthogonal projector matrices corresponding to all
strata of the experiment.
}
\details{
The strata decomposition is performed within this function. The first step
is to convert the list of block design matrices generated by
\code{\link{makeBlkDesMat}} to projection matrices using
\code{\link{projMat}}. The second step is to use these projection matrices
to project the raw data vector from one stratum to next stratum of the
experiment; the resulting matrix corresponds to each stratum is the
orthogonal projector matrix of the given stratum.
}
\examples{

design1 <- local({ 
  Ani = as.factor(LETTERS[c(1,2,3,4,
                            5,6,7,8)])
  Trt = as.factor(letters[c(1,1,1,1,
                            2,2,2,2)])
  data.frame(Ani, Trt, stringsAsFactors = TRUE )
})

blk.str = "Ani"
    
rT = terms(as.formula(paste("~", blk.str, sep = "")), keep.order = TRUE) 
blkTerm = attr(rT,"term.labels")
     
Z = makeBlkDesMat(design1, blkTerm)
Pb = makeOrthProjectors(Z)

}
\author{
Kevin Chang
}
