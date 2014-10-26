\name{toLatexTable}
\alias{toLatexTable}
\title{
Convert the R output to Latex Table
}
\description{
Print the Latex scripts on the screen for the user to output the table from the Latex output.
}
\usage{
toLatexTable(ANOVA, EF, fixed.names)
}
\arguments{
  \item{ANOVA}{
a matrix containing the coefficients of the variance components in EMS of ANOVA table generated by \code{\link{getCoefVC.onePhase}} or \code{\link{getCoefVC.twoPhase}}.
}
 \item{EF}{
a matrix containing the coefficient of the fixed effects components and the treatment average efficiency factors generated by \code{\link{getFixedEF.onePhase}} or \code{\link{getFixedEF.onePhase}} function.
}
  \item{fixed.names}{
a vector of character allows the users to modify symbols for the fixed effects.
}
}
\details{
Once the Latex script is generated, it requires the user to install and load two Latex packages: \code{booktabs} and \code{bm} to compile the Latex script.
}
\author{
Kevin Chang
}
\examples{
design1 <- local({ 
  Ani = as.factor(LETTERS[c(1,2,3,4,
                            5,6,7,8)])
  Trt = as.factor(letters[c(1,1,1,1,
                            2,2,2,2)])
  data.frame(Ani, Trt)
})

blk.str <- "Ani"
    
rT <- terms(as.formula(paste("~", blk.str, sep = "")), keep.order = TRUE) 
blkTerm <- attr(rT,"term.labels")
     
Z <- makeBlkDesMat(design1, blkTerm)


trt.str = "Trt"              
fT <- terms(as.formula(paste("~", trt.str, sep = "")), keep.order = TRUE)  #fixed terms

trtTerm <- attr(fT, "term.labels")
effectsMatrix <- attr(fT, "factor")        

T <- makeContrMat(design1, trtTerm, effectsMatrix, contr.vec = NA)

N <- makeOverDesMat(design1, trtTerm)

PNTginvATNP <- lapply( makeOrthProjectors(Z), function(z) infoDecompMat(z, T, N))

#Now construct variance matrices
Pb <- PNTginvATNP[sort(1:length(PNTginvATNP), decreasing=TRUE)]

v.mat <- getVMat.onePhase(Z.Phase1 = Z, design.df = design.df, var.comp = NA)
    
ANOVA <- getCoefVC.onePhase(Pb = Pb, design.df = design1, v.mat = v.mat, response = NA, 
	table.legend = FALSE, decimal = FALSE, digits = 2)
		
trt.Coef <- getTrtCoef(design1, trtTerm)

Replist = getTrtRep(design1, trtTerm)   
 
Rep <- Replist$Rep
trt.Sca <- Replist$Sca
    
effFactors = lapply(makeOrthProjectors(Z), function(z) 
      getEffFactor(z, T, N, Rep, trt.Sca))

effFactors <- effFactors[sort(1:length(effFactors), decreasing=TRUE)]

EF <- getFixedEF.onePhase(effFactors = effFactors, trt.Coef = trt.Coef,  T = T, Rep = Rep, 
	table.legend = FALSE, decimal = FALSE, digits = 2)

toLatexTable(ANOVA = ANOVA, EF = EF, fixed.names = c("\\\tau"))
}
