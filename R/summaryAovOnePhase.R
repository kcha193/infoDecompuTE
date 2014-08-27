#' Summarize an Theoretical Analysis of Variance Model of Single-Phase
#' Experiments
#' 
#' Computes the coefficients of the variance components for the expected mean
#' squares for single-phase experiments. The function accepts a data frame of
#' the experimetnal design with the structural formulae of the block and
#' treatment factors. Two tables containing the variance components of the
#' random effects and fixed effects are returned.
#' 
#' 
#' @param design.df a data frame containing the experimental design. Requires
#' every column be a \code{\link{factor}}.
#' @param blk.str a single string of characters containing the structural
#' formula for the block factors using the Wilkinson-Rogers' syntax.
#' @param trt.str a single string of characters containing the structural
#' formula for the treatment factors using the Wilkinson-Rogers' syntax.
#' @param var.comp a vector of characters containing the variance components of
#' interest this allows the user to specify the variance components to be shown
#' on the ANOVA table. This also allows the user to specify artificial stratum
#' to facilitate decomposition. Default is \code{NA}, which uses every random
#' factor as the variance components from \code{random.terms}.
#' @param trt.contr a list of treatment contrast vectors, this allows the user
#' to specify the contrasts for each treatment factor. Note that if this
#' arguement is used, it is necessary to specify the contrasts for every
#' treatment factor with the same order as \code{fixed.terms}. Default is
#' \code{NA}, which uses the C matrix described by John and Williams (1987).
#' @param table.legend a logical allows the users to use the legend for the
#' variance components of the ANOVA table for a large design. Default is
#' \code{FALSE}, which uses the original names.
#' @param response a numeric vector contains the responses from the experiment.
#' @param latex a logical allows the users to output the Latex script to Latex
#' table. Once the Latex script is generated, it requires the user to install
#' and load two Latex packages: \code{booktabs} and \code{bm} to complie the
#' Latex script.
#' @param fixed.names a vector of character allows the users to modify symbols
#' for the fixed effects for the Latex outputs.
#' @param decimal a logical allows users to display the coefficients as the
#' decimals. Default is \code{FALSE}, resulting in the use of \code{fractions}.
#' @param digits a integer indicating the number of decimal places. Default is
#' 2, resulting in 2 decimal places.
#' @return The values returned depends on the value of the \code{table.legend}
#' arguement. If \code{table.legend = FALSE}, this function will return a list
#' of two data frames. The first data frame contains the random effects and the
#' second data frame contains the fixed effects. If the \code{table.legend}
#' arguement is \code{TRUE}, then it will return a list containing two lists.
#' The first list consists of a data frame of random effects and a character
#' string for the legend. The second list consists of a data frame of fixed
#' effects and a character string for the legend.  If \code{response} arguement
#' is used, the random effect table will have one extra column with of mean
#' squares computed from the responses from the experiment.
#' @author Kevin Chang
#' @seealso \code{\link{terms}} for more information on the structural formula.
#' @references John J, Williams E (1987). \emph{Cyclic and computer generated
#' Designs}. Second edition. Chapman & Hall.
#' 
#' Nelder JA (1965b). "The Analysis of Randomized Experiments with Orthogonal
#' Block Structure. II. Treatment Structure and the General Analysis of
#' Variance." \emph{Proceedings of the Royal Society of London. Series A,
#' Mathematical and Physical Sciences}, 283(1393), 163-178.
#' 
#' Wilkinson GN, Rogers CE (1973). "Symbolic Description of Factorial Models
#' for Analysis of Variance." \emph{Applied Statistics}, 22(3), 392-399.
#' @keywords design
#' @examples
#' 
#' design1 <- local({ 
#'   Ani = as.factor(LETTERS[c(1,2,3,4,
#'                             5,6,7,8)])
#'   Trt = as.factor(letters[c(1,1,1,1,
#'                             2,2,2,2)])
#'   data.frame(Ani, Trt)
#' })
#' 
#' summaryAovOnePhase(design1, blk.str = "Ani", trt.str = "Trt") 
#' 
#' summaryAovOnePhase(design1, blk.str = "Ani", trt.str = "Trt", 
#' latex = TRUE, fixed.names = c("\\tau"))
#' 
#' 
#' @export summaryAovOnePhase
summaryAovOnePhase <- function(design.df, blk.str, trt.str, var.comp = NA, trt.contr = NA, 
    table.legend = FALSE, response = NA, latex = FALSE, fixed.names = NA, decimal = FALSE, 
    digits = 2) {
    
    design.df <- data.frame(lapply(design.df, factor))
    ######################################################################################### Main methods starts here-> Extract the fixed and random terms
    
    rT <- terms(as.formula(paste("~", blk.str, sep = "")), keep.order = TRUE)  #random terms
    
    rT.terms <- attr(rT, "term.labels")
    
    fT <- terms(as.formula(paste("~", trt.str, sep = "")), keep.order = TRUE)  #fixed terms
    
    ######################################################################################### Preparing the block structures
    #browser()
    Z <- makeBlkDesMat(design.df, rev(rT.terms))
    
    Pb <- makeOrthProjectors(Z)
    
	if(length(rT.terms) > 1)
		rT.terms = adjustEffectNames(effectsMatrix = attr(rT, "factors"), effectNames = rT.terms)
	
    if (names(Pb)[1] == "e") {
        names(Pb)[1] <- paste("Within", paste(unique(unlist(strsplit(names(Pb)[-1], 
            "[[:punct:]]"))), collapse = "."))
        names(Pb)[-1] <- rev(rT.terms)
        
    } else {
        names(Pb) <- rev(rT.terms)
        
    }
    names(Z)[-1] <- rev(rT.terms)
    
    ######################################################################################### Prepating the treatment structures
    
    trtTerm <- attr(fT, "term.labels")
    effectsMatrix <- attr(fT, "factor")
    #browser()
    
	if(length(trtTerm) > 1)
	  trtTerm = adjustEffectNames(effectsMatrix, trtTerm)
		
  	T <- makeContrMat(design.df = design.df, effectNames = trtTerm, effectsMatrix = effectsMatrix, 
  	                  contr.vec = trt.contr)
  	N <- makeOverDesMat(design.df = design.df, effectNames = trtTerm)
  	Rep <- getTrtRep(design.df, trtTerm)
  	trt.Coef <- getTrtCoef(design.df, trtTerm)
	
	#browser()
	
	#################################################################################################
	#Needs to check
	#allZero = apply(N, 2, function(x) all(x==0))
	#N = N[,!allZero]
	#T = lapply(T, function(x) x[!allZero, !allZero])
	#Rep = Rep[!allZero, ]
	
	
	################################################################################################
	
	#When there are treatment contrasts defined 
  	if (any(grepl("\\.", names(T)))) {
  	  colnames(Rep) <- trtTerm
  	  names(trt.Coef) <- trtTerm 
  	  Rep <- Rep[, sapply(strsplit(names(T), "\\."), function(x) x[1])]
  	  trt.Coef <- trt.Coef[sapply(strsplit(names(T), "\\."), function(x) x[1])]
  	} else {
  	 
  	  colnames(Rep) <- trtTerm
  	  names(trt.Coef) <- trtTerm      
	}
	
    ######################################################################################### Start calculating the VCs 1-phase experiment
    #browser()
    # pre- and post-multiply NTginvATN by block projection matrices
    PNTginvATNP <- lapply(Pb, function(z) infoDecompMat(z, T, N))
    
	
    # Now construct variance matrices
    PNTginvATNP <- PNTginvATNP[sort(1:length(PNTginvATNP), decreasing = TRUE)]
    
    # browser()
    v.mat <- getVMat.onePhase(Z.Phase1 = Z, design.df = design.df, var.comp = var.comp)
    
    if (all(is.na(var.comp))) {
        names(v.mat)[-1] <- rev(rT.terms)
    }
    # browser()
    
    ANOVA <- getCoefVC.onePhase(Pb = PNTginvATNP, design.df = design.df, v.mat = v.mat, 
        response = response, table.legend = table.legend, decimal = decimal, digits = digits)
    
    ############################################################################################################## 
    
    
    effFactors <- lapply(Pb, function(z) getEffFactor(z, T, N, Rep))
    effFactors <- effFactors[sort(1:length(effFactors), decreasing = TRUE)]
    
    EF <- getFixedEF.onePhase(effFactors = effFactors, trt.Coef = trt.Coef, T = T, Rep = Rep, 
        table.legend = table.legend, decimal = decimal, digits = digits)
    if (latex) {
        return(toLatexTable(ANOVA = ANOVA, EF = EF, fixed.names = fixed.names))
    } else {
        return(list(ANOVA = ANOVA, EF = EF))
    }
} 
