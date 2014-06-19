\name{projMat}
\alias{projMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct a Projection Matrix
}
\description{
Compute the projection matrix from a square matrix.
}
\usage{
projMat(X)
}
\arguments{
  \item{X}{
a square matrix.
}
}
\value{
A square matrix.
}

\author{
Kevin Chang
}

\examples{
m = matrix(1, nrow = 10, ncol = 3)
projMat(m) 
}
