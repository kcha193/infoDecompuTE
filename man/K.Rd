\name{K}
\alias{K}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Averaging Matrix 
}
\description{
Construct a n-by-n averaging matrix.}
\usage{
K(n)
}
\arguments{
  \item{n}{
a numeric describes the dimension of the averaging matrix.
}
}
\value{
This function returns a \eqn{n \times n} square matrix with all elements equal 1/n. 
}
\references{
John J, Williams E (1987). \emph{Cyclic and computer generated Designs}. Second edition. Chapman
& Hall.
}
\author{
Kevin Chang
}

\examples{
K(10)
}
