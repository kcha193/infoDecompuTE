summaryAovOnePhase <- function(design.df, blk.str, trt.str, var.comp = NA, trt.contr = NA, 
    table.legend = FALSE, response = NA, latex = FALSE, fixed.names = NA, decimal = FALSE, 
    digits = 2) {
    
    design.df <- data.frame(lapply(design.df, factor))
	
	newTerms = adjustMissingLevels(design.df, trt.str)
	#browser()	
	design.df = newTerms$design.df
	trt.str = newTerms$str.for
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
    
	effectsMatrix[nrow(effectsMatrix), effectsMatrix[nrow(effectsMatrix),]==2] <- 1
	
	if(length(trtTerm) > 1)
	  trtTerm = adjustEffectNames(effectsMatrix, trtTerm)
		
  	T <- makeContrMat(design.df = design.df, effectNames = trtTerm, effectsMatrix = effectsMatrix, 
  	                  contr.vec = trt.contr)
  	N <- makeOverDesMat(design.df = design.df, effectNames = trtTerm)
  	Replist <- getTrtRep(design.df, trtTerm)
	Rep <- Replist$Rep
  	trt.Sca <- Replist$Sca
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
     #browser()
    
    ANOVA <- getCoefVC.onePhase(Pb = PNTginvATNP, design.df = design.df, v.mat = v.mat, 
        response = response, table.legend = table.legend, decimal = decimal, digits = digits)
    
    ############################################################################################################## 
    
	#browser()
    effFactors <- lapply(Pb, function(z) getEffFactor(z, T, N, Rep, trt.Sca))
    effFactors <- effFactors[sort(1:length(effFactors), decreasing = TRUE)]
    
	#browser()
	
    EF <- getFixedEF.onePhase(effFactors = effFactors, trt.Coef = trt.Coef, T = T, Rep = Rep, 
        table.legend = table.legend, decimal = decimal, digits = digits)
    if (latex) {
        return(toLatexTable(ANOVA = ANOVA, EF = EF, fixed.names = fixed.names))
    } else {
        return(list(ANOVA = ANOVA, EF = EF))
    }
} 
