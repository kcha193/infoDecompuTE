summaryAovTwoPhase <- function(design.df, blk.str1, blk.str2, trt.str, var.comp = NA, blk.contr = NA, trt.contr = NA, 
    table.legend = FALSE, response = NA, latex = FALSE, fixed.names = NA, 
    decimal = FALSE, digits = 2){
    
    design.df = data.frame(lapply(design.df, factor))
    
    # Extract the fixed and random terms
    
    rT1 = terms(as.formula(paste("~", blk.str1, sep = "")), keep.order = TRUE)  #random terms phase 1
    rT2 = terms(as.formula(paste("~", blk.str2, sep = "")), keep.order = TRUE)  #random terms phase 2
    fT = terms(as.formula(paste("~", trt.str, sep = "")), keep.order = TRUE)  #fixed terms
    
    
    # Preparing the block
    # structures write('1.
    # Preparing the block
    # structure.', '')
    
    blkTerm1 = attr(rT1, "term.labels")
    blkTerm2 = attr(rT2, "term.labels")
    
    #check for complete confounding and changing the names of the block factors
    #browser()
    if(any(grepl("\\:", blkTerm2))){
          
      check.blkTerm2 = unique(unlist(strsplit(blkTerm2, "\\:")))
    } else {
      check.blkTerm2 = blkTerm2
    }
    
    if(any(grepl("\\:", blkTerm1))){
      check.blkTerm1 = unique(unlist(strsplit(blkTerm1, "\\:")))
    } else {
      check.blkTerm1 = blkTerm1
    }
   
   #browser()
    for(i in 1:length(check.blkTerm2)){
     if(length(check.blkTerm1) == 1){
         if(all(as.numeric(as.factor(design.df[, check.blkTerm1])) == as.numeric(as.factor(design.df[, check.blkTerm2[i]])))){
            cat("Note: Complete confounding between ", check.blkTerm1, " and ", check.blkTerm2[i], "!\n", sep = "")
            
           colnames(design.df)[which(colnames(design.df) == check.blkTerm1)] = 
                    paste(colnames(design.df)[which(colnames(design.df) == check.blkTerm1)], "CCW", sep = "")
           
           blk.str1 =  paste(blkTerm1[which(blkTerm1 == check.blkTerm1)], "CCW", sep = "")
           
           if(!is.na(blk.contr)){
              names(blk.contr)[which(names(blk.contr) == check.blkTerm1)] = 
                      paste(names(blk.contr)[which(names(blk.contr) == check.blkTerm1)], "CCW", sep = "")
           }
           check.blkTerm1 = paste(check.blkTerm1, "CCW", sep = "")
 
         }         
     }  else {
        check.temp = apply(design.df[, check.blkTerm1], 2, function(x) 
              all(as.numeric(as.factor(x))== as.numeric(as.factor(design.df[, check.blkTerm2[i]]))) )                                                      
          if(any(check.temp)){
             cat("Note: Complete confounding between ", check.blkTerm1[which(check.temp)], " and ", check.blkTerm2[i], "!\n", sep = "")              
            colnames(design.df)[which(colnames(design.df) == check.blkTerm1[which(check.temp)])] = 
                paste(colnames(design.df)[which(colnames(design.df) == check.blkTerm1[which(check.temp)])], 
                                                                                              "CCW", sep = "") 
            
            blk.str1 =  gsub(check.blkTerm1[which(check.temp)],  paste(check.blkTerm1[which(check.temp)], "CCW", sep = ""), blk.str1)      
            if(!is.na(blk.contr)){
              names(blk.contr)[which(names(blk.contr) == check.blkTerm1[which(check.temp)])] = 
                      paste(names(blk.contr)[which(names(blk.contr) == check.blkTerm1[which(check.temp)])], 
                            "CCW", sep = "")
           }      
            check.blkTerm1[which(check.blkTerm1 == check.blkTerm1[which(check.temp)])] = 
                  paste(check.blkTerm1[which(check.blkTerm1 == check.blkTerm1[which(check.temp)])], "CCW", sep = "")
          }
      }
    }
    
    
    rT1 = terms(as.formula(paste("~", blk.str1, sep = "")), keep.order = TRUE)  #random terms phase 1
    blkTerm1 = attr(rT1, "term.labels")

    
    Z1 = makeBlkDesMat(design.df, rev(blkTerm1))
    Z2 = makeBlkDesMat(design.df, rev(blkTerm2))
    
    # write('2. Defining the block structures of second Phase.', '')
    Pb <- makeOrthProjectors(Z2)
    
  
checkBrack = function(x, str.for) {
        #browser()
          str.for = gsub("/", "\\:",str.for)
         #storeBrack = unlist(strsplit(gsub("([[:punct:]])", "\\.\\1\\.",  str.for), "\\."))
         storeBrack = unlist(strsplit(unlist(strsplit(str.for, "\\(")), "\\)"))
          
        if(any(storeBrack==""))
          storeBrack = storeBrack[-which(storeBrack=="")]
          
          if(any(x == storeBrack)) return(x)         
          
          nameX =  unlist(strsplit(x, "([[:punct:]])"))
          
          x = unlist(strsplit(gsub("([[:punct:]])", "\\.\\1\\.", x), "\\."))
          

          combName = ""
          openBrack = FALSE
          openBrack1 = FALSE
          while(length(x) > 0){
            if(grepl(x[1], "[[:punct:]]")){
              combName = paste(combName, x[1], sep = "") 
            } else{
               matchBrack = storeBrack[grep(x[1],storeBrack)]
            
               eachMatchBrack = unlist(strsplit(matchBrack, "[[:punct:]]"))
               
               if(any(eachMatchBrack==""))
                eachMatchBrack = eachMatchBrack[-which(eachMatchBrack=="")]
          
               if(length(eachMatchBrack) == 1){
                if(grepl("\\:", matchBrack) &&
                 any(apply(sapply(nameX[-1], function(x) grepl(x, storeBrack)), 1, all))){
                  combName = paste(combName, x[1], sep = "") 
               
                }else if(grepl("\\:", matchBrack) &&
                    grepl("\\*", storeBrack[grep(nameX[-1][1],storeBrack)])){
                      
                    combName = paste(combName, "(", x[1], sep = "") 
                    openBrack = TRUE
                } else {
                               
                combName = paste(combName, x[1], sep = "")
                } 
               } else {
                  if((which(eachMatchBrack==x[1]) == 1) && 
                              any(!is.na(match(eachMatchBrack, x[-1])))){
                     combName = paste(combName, "(", x[1], sep = "") 
                      openBrack1 = TRUE
                   } else if((which(eachMatchBrack==x[1]) == length(eachMatchBrack)) || 
                                  all(is.na(match(eachMatchBrack, x[-1])))){ 
                     combName = paste(combName, x[1], ")", sep = "") 
                     
                     if(openBrack && openBrack1){
                        combName = paste(combName, ")", sep = "")
                        openBrack = FALSE 
                         openBrack1 = FALSE
                     }
                     
                  }else {
                     combName = paste(combName, x[1], sep = "") 
                  }
               }
            }
            x = x[-1]
          }
          
          return(combName)
     }

     checkCross = attr(rT2, "factor")
     
    if (any(grepl(":", blkTerm2))) {
        check.blkTerm2 = blkTerm2[grep(":", blkTerm2)]
        
        checkCross = checkCross[, check.blkTerm2]
        
        if (!is.matrix(checkCross)) {
            temp = matrix(checkCross)
            rownames(temp) = names(checkCross)
            colnames(temp) = check.blkTerm2
            checkCross = temp
        }
        
        split.check.blkTerm2 = strsplit(check.blkTerm2, ":")
        
        for (i in 1:length(split.check.blkTerm2)) {
            newName = ""
            for (j in 1:(length(split.check.blkTerm2[[i]]) - 1)) {
                if (checkCross[split.check.blkTerm2[[i]][j], i] == 1) {
                  newName = paste(newName, split.check.blkTerm2[[i]][j], "*", sep = "")
                } else {
                  newName = paste(newName, split.check.blkTerm2[[i]][j], ":", sep = "")
                }
            }
            newName = paste(newName, split.check.blkTerm2[[i]][length(split.check.blkTerm2[[i]])], sep = "")
            check.blkTerm2[i] = newName
        }
         #browser()
        
         len = intersect(grep("\\:", check.blkTerm2), grep("\\*", check.blkTerm2))
         if(length(len)> 0){
         for(i in 1:length(len))
            check.blkTerm2[len[i]] = 
                  checkBrack(check.blkTerm2[len[i]], blk.str2)
                  
        }
                  
        blkTerm2[grep(":", blkTerm2)] = check.blkTerm2
    }
    
    
    if (names(Pb)[1] == "e") {
        names(Pb)[1] = paste("Within", paste(unique(unlist(strsplit(names(Pb)[-1], "[[:punct:]]+"))), collapse = "."))
        
        names(Pb)[-1] = rev(blkTerm2)
    } else {
        names(Pb) = rev(blkTerm2)
    }     
    

    # write('3. Defining the block structures of first phase within second Phase.', '')
    blkTerm1 = attr(rT1, "term.labels")
    effectsMatrix = attr(rT1, "factor")
             
  checkCross = attr(rT1, "factor")
    if (any(grepl(":", blkTerm1))) {
        check.blkTerm1 = blkTerm1[grep(":", blkTerm1)]
        
        checkCross = checkCross[, check.blkTerm1]
        
        if (!is.matrix(checkCross)) {
            temp = matrix(checkCross)
            rownames(temp) = names(checkCross)
            colnames(temp) = check.blkTerm1
            checkCross = temp
        }
        
        split.check.blkTerm1 = strsplit(check.blkTerm1, ":")
        
        for (i in 1:length(split.check.blkTerm1)) {
            newName = ""
            for (j in 1:(length(split.check.blkTerm1[[i]]) - 1)) {
                if (checkCross[split.check.blkTerm1[[i]][j], i] == 1) {
                  newName = paste(newName, split.check.blkTerm1[[i]][j], "*", sep = "")
                } else {
                  newName = paste(newName, split.check.blkTerm1[[i]][j], ":", sep = "")
                }
            }
            newName = paste(newName, split.check.blkTerm1[[i]][length(split.check.blkTerm1[[i]])], sep = "")
            check.blkTerm1[i] = newName
        }
 
         len = intersect(grep("\\:", check.blkTerm1), grep("\\*", check.blkTerm1))
         if(length(len)> 0){
         for(i in 1:length(len))
            check.blkTerm1[len[i]] = 
                  checkBrack(check.blkTerm1[len[i]], blk.str1)
                  
        }
  
        blkTerm1[grep(":", blkTerm1)] =check.blkTerm1
        colnames(effectsMatrix) = blkTerm1
    }
 
		 
    T = makeContrMat(design.df, effectNames = blkTerm1, effectsMatrix = effectsMatrix, contr.vec = blk.contr)
    N = makeOverDesMat(design.df, blkTerm1)
    
     
    #browser()
     
    res = paste("Within", paste(unique(unlist(strsplit(names(T), "[[:punct:]]+"))), collapse = "."))
    
    Pb1 <- lapply(Pb, function(z) infoDecompMat(z, T, N))
    
    # t.name = unique(unlist(strsplit(names(T), '[[:punct:]]')))
    t.name = names(T)
    
	if(length(Pb) > 1){
		pb1.names = lapply((Pb1[-1]), names)
    }else{
		pb1.names = lapply(Pb1, names)
	}
	
    # browser()
    for (i in 1:length(pb1.names)) {
        comp = t.name %in% pb1.names[[i]]
        if (any(comp)) {
            break
        } else if (i == length(pb1.names)) {
            names(Pb1)[1] = ""
        }
    }
    
    for (i in 1:length(Pb1)) {
        names(Pb1[[i]])[which(names(Pb1[[i]]) == "Residual")] = res
    }
    
    # Preparing the treatment
    # structures write('4.
    # Preparing the treatment
    # structure.', '')
    
    trtTerm = attr(fT, "term.labels")
    effectsMatrix = attr(fT, "factor")
    
   
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
    
    colnames(Rep) = names(T)
    names(trt.Coef) = names(T)
    
    # Start calculating the
    # VCs 2-phase experiment
    # write('5. Start
    # calculating the
    # variance components.',
    # '') write('6. Pre- and
    # post-multiply NTginvATN
    # by block projection
    # matrices.', '')
    
    
    PNTginvATNP <- lapply(Pb1, function(y) lapply(y, function(z) infoDecompMat(z, T, N)))
    
    PNTginvATNP <- PNTginvATNP[sort(1:length(PNTginvATNP), decreasing = TRUE)]
    
    
    v.mat = getVMat.twoPhase(Z.Phase1 = Z1, Z.Phase2 = Z2, design.df = design.df, var.comp = var.comp)
    
    if (all(is.na(var.comp))) {
        names(v.mat)[-1] = c(rev(blkTerm1), rev(blkTerm2))
    }
    
    
    ANOVA = getCoefVC.twoPhase(Pb = PNTginvATNP, design.df = design.df, v.mat = v.mat, response, table.legend, decimal= decimal, digits = digits)
    
    
    effFactors = lapply(Pb1, function(y) lapply(y, function(z) getEffFactor(z, T, N, Rep)))
    effFactors <- effFactors[sort(1:length(effFactors), decreasing = TRUE)]
    
    
    EF = getFixedEF.twoPhase(effFactors = effFactors,  trt.Coef = trt.Coef, T = T, Rep = Rep, table.legend, decimal= decimal, digits = digits)
    
    if (latex) {
        return(toLatexTable(ANOVA = ANOVA, EF = EF, fixed.names = fixed.names))
    } else {
        return(list(ANOVA = ANOVA, EF = EF))
    }
    
    
} 
