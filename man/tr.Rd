\name{tr}
\alias{tr}
\title{
Trace of the Matrix
}
\description{
Compute the trace of the square matrix.
}
\usage{
tr(X)
}
\arguments{
  \item{X}{
a square matrix.
}
}

\value{
A numeric value.
}
\references{
John J, Williams E (1987). \emph{Cyclic and computer generated Designs}. Second edition. Chapman
& Hall.
}
\author{
Kevin
}

\seealso{
 \code{\link{diag}}
}
\examples{
m = matrix(1, nrow = 10, ncol = 10)
tr(m)   
}
