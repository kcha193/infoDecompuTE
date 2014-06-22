summaryAovOnePhase <- 
function(design.df, blk.str, trt.str, var.comp = NA, trt.contr = NA, table.legend = FALSE, 
	response = NA, latex = FALSE, fixed.names = NA, decimal = FALSE, digits = 2) {
    
    design.df = data.frame(lapply(design.df, factor))
    ######################################################################################### Main methods starts here-> Extract the fixed and random terms
    
    rT = terms(as.formula(paste("~", blk.str, sep = "")), keep.order = TRUE)  #random terms
    
    rT.terms = attr(rT, "term.labels")
    
    fT = terms(as.formula(paste("~", trt.str, sep = "")), keep.order = TRUE)  #fixed terms
    
    ######################################################################################### Preparing the block structures
    
    Z = makeBlkDesMat(design.df, rev(rT.terms))
    
    Pb <- makeOrthProjectors(Z)
        
    checkBrack = function(x, str.for) {
        # browser()
        str.for = gsub("/", "\\:", str.for)
        # storeBrack = unlist(strsplit(gsub('([[:punct:]])', '\\.\\1\\.', str.for), '\\.'))
        storeBrack = unlist(strsplit(unlist(strsplit(str.for, "\\(")), "\\)"))
        
        if (any(storeBrack == "")) 
            storeBrack = storeBrack[-which(storeBrack == "")]
        
        if (any(x == storeBrack)) 
            return(x)
        
        nameX = unlist(strsplit(x, "([[:punct:]])"))
        
        x = unlist(strsplit(gsub("([[:punct:]])", "\\.\\1\\.", x), "\\."))
        
        
        combName = ""
        openBrack = FALSE
        openBrack1 = FALSE
        while (length(x) > 0) {
            if (grepl(x[1], "[[:punct:]]")) {
                combName = paste(combName, x[1], sep = "")
            } else {
                matchBrack = storeBrack[grep(x[1], storeBrack)]
                eachMatchBrack = unlist(strsplit(matchBrack, "[[:punct:]]"))
                
                if (any(eachMatchBrack == "")) 
                  eachMatchBrack = eachMatchBrack[-which(eachMatchBrack == "")]
                
                if (length(eachMatchBrack) == 1) {
                  if (grepl("\\:", matchBrack) && any(apply(sapply(nameX[-1], function(x) grepl(x, storeBrack)), 1, 
                    all))) {
                    combName = paste(combName, x[1], sep = "")
                    
                  } else if (grepl("\\:", matchBrack) && grepl("\\*", storeBrack[grep(nameX[-1][1], storeBrack)])) {
                    
                    combName = paste(combName, "(", x[1], sep = "")
                    openBrack = TRUE
                  } else {
                    
                    combName = paste(combName, x[1], sep = "")
                  }
                } else {
                  if ((which(eachMatchBrack == x[1]) == 1) && any(!is.na(match(eachMatchBrack, x[-1])))) {
                    combName = paste(combName, "(", x[1], sep = "")
                    openBrack1 = TRUE
                  } else if ((which(eachMatchBrack == x[1]) == length(eachMatchBrack)) || all(is.na(match(eachMatchBrack, 
                    x[-1])))) {
                    combName = paste(combName, x[1], ")", sep = "")
                    
                    if (openBrack && openBrack1) {
                      combName = paste(combName, ")", sep = "")
                      openBrack = FALSE
                      openBrack1 = FALSE
                    }
                    
                  } else {
                    combName = paste(combName, x[1], sep = "")
                  }
                }
            }
            x = x[-1]
        }
        
        return(combName)
    }
    
    
    checkCross = attr(rT, "factor")
    
    if (any(grepl(":", rT.terms))) {
        check.rT.terms = rT.terms[grep(":", rT.terms)]
        
        checkCross = checkCross[, check.rT.terms]
        
        if (!is.matrix(checkCross)) {
            temp = matrix(checkCross)
            rownames(temp) = names(checkCross)
            colnames(temp) = check.rT.terms
            checkCross = temp
        }
        
        split.check.rT.terms = strsplit(check.rT.terms, ":")
        
        for (i in 1:length(split.check.rT.terms)) {
            newName = ""
            for (j in 1:(length(split.check.rT.terms[[i]]) - 1)) {
                if (checkCross[split.check.rT.terms[[i]][j], i] == 1) {
                  newName = paste(newName, split.check.rT.terms[[i]][j], "*", sep = "")
                } else {
                  newName = paste(newName, split.check.rT.terms[[i]][j], ":", sep = "")
                }
            }
            newName = paste(newName, split.check.rT.terms[[i]][length(split.check.rT.terms[[i]])], sep = "")
            check.rT.terms[i] = newName
        }
        
        len = intersect(grep("\\:", check.rT.terms), grep("\\*", check.rT.terms))
        if (length(len) > 0) {
            for (i in 1:length(len)) check.rT.terms[len[i]] = checkBrack(check.rT.terms[len[i]], blk.str)
            
        }
        
        
        rT.terms[grep(":", rT.terms)] = check.rT.terms
    }
    # browser()
    
    if (names(Pb)[1] == "e") {
        names(Pb)[1] = paste("Within", paste(unique(unlist(strsplit(names(Pb)[-1], "[[:punct:]]+"))), collapse = "."))
        names(Pb)[-1] = rev(rT.terms)
        
    } else {
        names(Pb) = rev(rT.terms)
        
    }
    names(Z)[-1] = rev(rT.terms)
    
    ######################################################################################### Prepating the treatment structures
    
    trtTerm = attr(fT, "term.labels")
    effectsMatrix = attr(fT, "factor")
    # browser()

	 if (any(grepl(":", trtTerm))) {
        check.trtTerm = trtTerm[grep(":", trtTerm)]
        
        effectsMatrix1 = effectsMatrix[, check.trtTerm]
        
        if (!is.matrix(effectsMatrix1)) {
            temp = matrix(effectsMatrix1)
            rownames(temp) = names(effectsMatrix1)
            colnames(temp) = check.trtTerm
            effectsMatrix1 = temp
        }
        
        split.check.trtTerm = strsplit(check.trtTerm, ":")
        
        for (i in 1:length(split.check.trtTerm)) {
            newName = ""
            for (j in 1:(length(split.check.trtTerm[[i]]) - 1)) {
                if (effectsMatrix1[split.check.trtTerm[[i]][j], i] == 1) {
                  newName = paste(newName, split.check.trtTerm[[i]][j], "*", sep = "")
                } else {
                  newName = paste(newName, split.check.trtTerm[[i]][j], ":", sep = "")
                }
            }
            newName = paste(newName, split.check.trtTerm[[i]][length(split.check.trtTerm[[i]])], sep = "")
            check.trtTerm[i] = newName
        }
        
        len = intersect(grep("\\:", check.trtTerm), grep("\\*", check.trtTerm))
        if (length(len) > 0) {
            for (i in 1:length(len)) check.trtTerm[len[i]] = checkBrack(check.trtTerm[len[i]], trt.str)
            
        }
        
        trtTerm[grep(":", trtTerm)] = check.trtTerm
        colnames(effectsMatrix) = trtTerm
        
    }
    
    T = makeContrMat(design.df = design.df, effectNames = trtTerm, effectsMatrix = effectsMatrix, 
		contr.vec = trt.contr)
    N = makeOverDesMat(design.df = design.df,  effectNames = trtTerm)
    Rep = getTrtRep(design.df, trtTerm)
    trt.Coef = getTrtCoef(design.df, trtTerm)
    
    if (any(grepl("\\.", names(T)))) {
        colnames(Rep) = trtTerm
        names(trt.Coef) = trtTerm
        
        Rep = Rep[, sapply(strsplit(names(T), "\\."), function(x) x[1])]
        trt.Coef = trt.Coef[sapply(strsplit(names(T), "\\."), function(x) x[1])]
    }
    # browser()
    
  
    colnames(Rep) = names(T)
    names(trt.Coef) = names(T)
    ######################################################################################### Start calculating the VCs 1-phase experiment
    
    # pre- and post-multiply NTginvATN by block projection matrices
    PNTginvATNP <- lapply(Pb, function(z) infoDecompMat(z, T, N))
    
    # Now construct variance matrices
    PNTginvATNP <- PNTginvATNP[sort(1:length(PNTginvATNP), decreasing = TRUE)]
    
    # browser()
    v.mat = getVMat.onePhase(Z.Phase1 = Z, design.df = design.df, var.comp = var.comp)
    
    if (all(is.na(var.comp))) {
        names(v.mat)[-1] = rev(rT.terms)
    }
    # browser()
    
    ANOVA = getCoefVC.onePhase(Pb = PNTginvATNP, design.df = design.df, v.mat = v.mat, response = response, 
			table.legend = table.legend, decimal= decimal, digits = digits)
    
    ############################################################################################################## 
    
    
    effFactors = lapply(Pb, function(z) getEffFactor(z, T, N, Rep))
    effFactors <- effFactors[sort(1:length(effFactors), decreasing = TRUE)]
    
    EF = getFixedEF.onePhase(effFactors = effFactors, trt.Coef = trt.Coef, T = T, Rep = Rep, table.legend = table.legend, 
		decimal= decimal, digits = digits)
    if (latex) {
        return(toLatexTable(ANOVA = ANOVA, EF = EF, fixed.names = fixed.names))
    } else {
        return(list(ANOVA = ANOVA, EF = EF))
    }
} 
