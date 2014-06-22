makeContrMat <- function(design.df, effectNames, effectsMatrix, contr.vec) {
    # Square contrast matrix with the specific contrast (trtContr)
    indMatrix1 <- function(x, n, trtContr) {
        if (x == 1) 
            X <- trtContr else if (x == 2) 
            X <- identityMat(n) else X <- K(n)
        return(X)
    }
    
    # Square contrast matrix without defining the specific contrast
    indMatrix <- function(x, n) {
        if (x == 1) 
            X <- J(n) else if (x == 2) 
            X <- identityMat(n) else X <- K(n)
        return(X)
    }
    
    # transform each treatment contrast matrix to C matrix
    transContrToT <- function(fact, contr) {
        T <- matrix(0, nrow = nlevels(fact), ncol = nlevels(fact))
        if (is.matrix(contr)) {
            rownames(contr) <- fact
            for (i in 1:ncol(contr)) {
                x <- contr[levels(fact), i]
                T <- T + (x %*% t(x))/as.numeric(t(x) %*% x)
            }
            
        } else {
            names(contr) <- fact
            x <- contr[levels(fact)]
            T <- T + (x %*% t(x))/as.numeric(t(x) %*% x)
        }
        return(T)
    }
    
    # obtaining the numbers of the levels from the design
    if (any(grepl("[[:punct:]]", effectNames))) {
        uniqueTrtCols <- unique(unlist(strsplit(effectNames, "[[:punct:]]")))
        nLevels <- sapply(design.df[, uniqueTrtCols], function(x) nlevels(as.factor(x)))
        
    } else if (length(effectNames) == 1) {
        nLevels <- nlevels(design.df[, effectNames])
        names(nLevels) <- effectNames
        
    } else {
        nLevels <- sapply(design.df[, effectNames], function(x) nlevels(as.factor(x)))
    }
    
    nLevels <- nLevels[rownames(effectsMatrix)]
    nEffects <- ncol(effectsMatrix)
    
    
    # Without the contrasts specifically defined
    if (all(is.na(contr.vec))) {
        X <- as.list(rep(1, nEffects))
        names(X) <- colnames(effectsMatrix)
        
        for (i in 1:nrow(effectsMatrix)) {
            matList <- lapply(effectsMatrix[i, ], function(y) indMatrix(y, nLevels[i]))
            for (j in 1:nEffects) X[[j]] <- X[[j]] %x% matList[[j]]
        }
        
    } else {
        new.contr.vec <- vector(length = length(effectNames), mode = "list")
        names(new.contr.vec) <- effectNames
        
        for (i in 1:length(new.contr.vec)) {
            new.contr.vec[i] <- ifelse(is.null(contr.vec[[names(new.contr.vec)[i]]]), 
                NA, contr.vec[names(new.contr.vec)[i]])
        }
        
        
        
        # Construct the interactions
        if (any(grepl("[[:punct:]]", effectNames))) {
            interTerms <- grep("[[:punct:]]", effectNames)
            
            for (i in interTerms) {
                mainTerms <- unlist(strsplit(effectNames[[i]], "[[:punct:]]"))
                mainTerms.list <- vector(length = length(mainTerms), mode = "list")
                names(mainTerms.list) <- mainTerms
                
                for (j in 1:length(mainTerms)) {
                  if (is.null(names(contr.vec[[mainTerms[j]]]))) 
                    mainTerms.list[[j]] <- 0 else mainTerms.list[[j]] <- names(contr.vec[[mainTerms[j]]])
                }
                
                
                new.contr.vec[[i]] <- vector(length = nlevels(interaction(mainTerms.list)), 
                  mode = "list")
                # names(contr.vec)[[i]] = effectNames[i]
                
                names(new.contr.vec[[i]]) <- levels(interaction(mainTerms.list))
            }
        }
        
        X <- list()
        count <- 1
        
        contr.vec <- new.contr.vec
        
        totalLength <- length(contr.vec[-which(sapply(contr.vec, class) == "list")]) + 
            sum(sapply(contr.vec[which(sapply(contr.vec, class) == "list")], length))
        
        newEffectsMatrix <- matrix(0, ncol = totalLength, nrow = nrow(effectsMatrix))
        rownames(newEffectsMatrix) <- rownames(effectsMatrix)
        
        # construct a new effects matrix to use for the Kronecker products
        for (i in 1:length(contr.vec)) {
            
            
            if (is.list(contr.vec[[i]])) {
                tmpX <- rep(1, length(contr.vec[[i]]))
                names(tmpX) <- paste(names(contr.vec[i]), names(contr.vec[[i]]), sep = ".")
                
                tmpCounter <- count:(count + length(contr.vec[[i]]) - 1)
                
                X <- c(X, tmpX)
                rowNames <- unlist(strsplit(names(contr.vec[i]), "[[:punct:]]"))
                newEffectsMatrix[rowNames, tmpCounter] <- 1
                count <- count + length(contr.vec[[i]])
                
            } else {
                tmpX <- 1
                names(tmpX) <- names(contr.vec[i])
                X <- c(X, tmpX)
                
                newEffectsMatrix[, count] <- effectsMatrix[, i]
                
                tmpCounter <- count <- count + 1
            }
        }
        
        colnames(newEffectsMatrix) <- names(X)
        
        for (i in 1:nrow(newEffectsMatrix)) {
            if (is.list(contr.vec[[rownames(effectsMatrix)[i]]])) {
                trtContr <- lapply(contr.vec[[rownames(effectsMatrix)[i]]], function(x) transContrToT(design.df[, 
                  rownames(effectsMatrix)[i]], x))
                
                matList <- vector(mode = "list", length = totalLength)
                
                for (j in 1:ncol(newEffectsMatrix)) {
                  
                  if (newEffectsMatrix[i, j] == 1) {
                    trtNames <- match(unlist(strsplit(colnames(newEffectsMatrix)[j], 
                      "\\.")), names(trtContr))
                    
                    matList[[j]] <- trtContr[[trtNames[-which(is.na(trtNames))]]]
                  } else if (newEffectsMatrix[i, j] == 2) {
                    matList[[j]] <- identityMat(nLevels[i])
                  } else {
                    matList[[j]] <- K(nLevels[i])
                  }
                }
            } else {
                if (all(is.na(contr.vec[[i]]))) {
                  matList <- lapply(newEffectsMatrix[i, ], function(y) indMatrix(y, 
                    nLevels[i]))
                  
                } else {
                  trtContr <- transContrToT(design.df[, rownames(effectsMatrix)[i]], 
                    contr.vec[[rownames(effectsMatrix)[i]]])
                  matList <- lapply(newEffectsMatrix[i, ], function(y) indMatrix1(y, 
                    nLevels[i], trtContr))
                  
                }
            }
            
            
            for (j in 1:length(X)) X[[j]] <- X[[j]] %x% matList[[j]]
            
        }
        
        
        
    }
    
    names(X) <- gsub("\\.0", "", names(X))
    return(X)
} 
