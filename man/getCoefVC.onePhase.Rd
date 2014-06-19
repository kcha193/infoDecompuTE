\name{getCoefVC.onePhase}
\alias{getCoefVC.onePhase}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCoefVC.onePhase(Pb, design.df, v.mat, response, table.legend, decimal, digits)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pb}{
%%     ~~Describe \code{Pb} here~~
}
  \item{design.df}{
%%     ~~Describe \code{design.df} here~~
}
  \item{v.mat}{
%%     ~~Describe \code{v.mat} here~~
}
  \item{response}{
%%     ~~Describe \code{response} here~~
}
  \item{table.legend}{
%%     ~~Describe \code{table.legend} here~~
}
  \item{decimal}{
%%     ~~Describe \code{decimal} here~~
}
  \item{digits}{
%%     ~~Describe \code{digits} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Pb, design.df, v.mat, response, table.legend, decimal, 
    digits) 
{
    if (all(is.na(response))) {
        response = rep(NA, nrow(design.df))
    }
    MS = lapply(Pb, function(q) lapply(q, function(w) t(response) \%*\% 
        w \%*\% response))
    V = v.mat
    VC <- rep("1", length(V) + 2)
    names(VC) = c("DF", names(V), "MS")
    VC = t(data.frame(VC))
    for (i in 1:length(Pb)) {
        tmp <- matrix(0, nrow = length(names(Pb[[i]])), ncol = length(V) + 
            1, dimnames = list(names(Pb[[i]]), c(names(V), "MS")))
        for (j in 1:(length(names(Pb[[i]])))) {
            for (z in 1:(length(V))) {
                tmp[j, z] <- tr(Pb[[i]][[j]] \%*\% V[[z]])
            }
            tmp[j, (length(V) + 1)] = as.numeric(MS[[i]][[j]])
        }
        if (nrow(tmp) == 1 && rownames(tmp) == "Residual") {
            tmp = c(tmp[1], tmp/tmp[1])
            if (decimal) {
                VC = rbind(VC, c(round(tmp[-length(tmp)], digits = digits), 
                  round(tmp[length(tmp)], digits = 3)))
            }
            else {
                VC = rbind(VC, c(attr(fractions(tmp[-length(tmp)]), 
                  "fracs"), round(tmp[length(tmp)], digits = 3)))
            }
            if (grepl("Within", names(Pb[i]))) {
                rownames(VC)[nrow(VC)] = paste(names(Pb[i]), 
                  sep = " ")
            }
            else {
                rownames(VC)[nrow(VC)] = paste("Between", names(Pb[i]), 
                  sep = " ")
            }
        }
        else {
            VC = rbind(VC, character(length = length(V) + 2))
            if (grepl("Within", names(Pb[i]))) {
                rownames(VC)[nrow(VC)] = paste(names(Pb[i]), 
                  sep = " ")
            }
            else {
                rownames(VC)[nrow(VC)] = paste("Between", names(Pb[i]), 
                  sep = " ")
            }
            rownames(tmp) = paste("  ", rownames(tmp), sep = " ")
            if (decimal) {
                tmp = t(apply(tmp, 1, function(x) c(round(c(x[1], 
                  x[-length(x)]/x[1]), digits = digits), as.character(round(x[length(x)]/x[1], 
                  digits = 3)))))
            }
            else {
                tmp = t(apply(tmp, 1, function(x) c(attr(fractions(c(x[1], 
                  x[-length(x)]/x[1])), "fracs"), as.character(round(x[length(x)]/x[1], 
                  digits = 3)))))
            }
            VC = rbind(VC, tmp)
        }
    }
    if (all(is.na(response))) {
        VC = VC[, -ncol(VC)]
    }
    if (length((which(apply(apply(VC[, -1], 2, function(x) x == 
        VC[, "e"]), 2, all)))) > 1) {
        VC = noquote(VC[-1, -2])
    }
    else {
        VC = noquote(VC[-1, ])
    }
    if (table.legend) {
        Legend = paste(paste(letters[1:(length(colnames(VC)) - 
            1)], colnames(VC)[-1], sep = " = "))
        colnames(VC)[-1] = letters[1:(length(colnames(VC)) - 
            1)]
        VC = list(VC = VC, Legend = Legend)
    }
    return(VC)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
