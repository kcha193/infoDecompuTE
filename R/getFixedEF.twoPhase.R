getFixedEF.twoPhase <- function(effFactors, trt.Coef, T, Rep, table.legend, decimal, 
    digits) {
    
    trt <- numeric(length(trt.Coef) + ncol(Rep))
    names(trt) <- c(names(T), paste("eff", names(T), sep = "."))
    noindent <- FALSE
    # browser()
    
    for (i in 1:length(effFactors)) {
        trt <- rbind(trt, character(length = length(T) * 2))
        if (grepl("Within", names(effFactors[i]))) {
            rownames(trt)[nrow(trt)] <- paste(names(effFactors[i]), sep = " ")
        } else if (names(effFactors[i]) == "") {
            trt <- trt[-(nrow(trt)), ]
            noindent <- TRUE
        } else {
            rownames(trt)[nrow(trt)] <- paste("Between", names(effFactors[i]), sep = " ")
        }
        
        for (k in 1:length(effFactors[[i]])) {
            trt <- rbind(trt, character(length = length(T) * 2))
            if (grepl("Within", names(effFactors[[i]][k])) && length(effFactors[[i]]) == 
                1) {
                trt <- trt[-(nrow(trt)), ]
                noindent <- TRUE
            } else if (grepl("Within", names(effFactors[[i]][k]))) {
                if (noindent) {
                  rownames(trt)[nrow(trt)] <- names(effFactors[[i]][k])
                } else {
                  rownames(trt)[nrow(trt)] <- paste("  ", names(effFactors[[i]][k]), 
                    sep = " ")
                }
            } else {
                if (noindent) {
                  rownames(trt)[nrow(trt)] <- paste("Between", names(effFactors[[i]][k]), 
                    sep = " ")
                } else {
                  rownames(trt)[nrow(trt)] <- paste("   Between", names(effFactors[[i]][k]), 
                    sep = " ")
                }
            }
            
            for (j in 1:length(effFactors[[i]][[k]])) {
                if (is.null(effFactors[[i]][[k]][[j]])) 
                  next
                
                if (decimal) {
                  
                  trt.temp <- round(c(trt.Coef * unlist(effFactors[[i]][[k]][[j]]), 
                    unlist(effFactors[[i]][[k]][[j]])), digits = digits)
                } else {
                  trt.temp <- attr(fractions(c(trt.Coef * unlist(effFactors[[i]][[k]][[j]]), 
                    unlist(effFactors[[i]][[k]][[j]]))), "fracs")
                }
                
                trt <- rbind(trt, trt.temp)
                
                if (noindent) {
                  rownames(trt)[nrow(trt)] <- paste("  ", names(effFactors[[i]][[k]][j]), 
                    sep = " ")
                } else {
                  rownames(trt)[nrow(trt)] <- paste("     ", names(effFactors[[i]][[k]][j]), 
                    sep = " ")
                }
            }
            
        }
        noindent <- FALSE
    }
    
    
    trt <- trt[-1, ]
    
    trt <- noquote(ifelse(trt == "NaN", "", trt))
    
    if (table.legend) {
        Legend <- paste(paste(letters[1:(length(colnames(trt)))], colnames(trt), sep = " = "))
        colnames(trt) <- letters[1:(length(colnames(trt)))]
        trt <- list(trt = trt, Legend = Legend)
    }
    
    
    return(trt)
    
} 
