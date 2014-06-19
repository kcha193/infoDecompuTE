\name{identityMat}
\alias{identityMat}
\title{
Identity Matrix
}
\description{
Construct an identity matrix.
}
\usage{
identityMat(n)
}
\arguments{
  \item{n}{
a numeric describes the dimension of the identity matrix.
}
}
\value{
This function returns a matrix with the diagonal elements equal to one and the off-diagonal elements equal to zero.
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
identityMat(10)
}
