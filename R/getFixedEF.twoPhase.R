getFixedEF.twoPhase <- function(effFactors, trt.Sca, T, Rep, table.legend, decimal, 
    digits) {
    
    trt <- numeric(length(trt.Sca) + ncol(Rep))
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
                
				effCoefList = effFactors[[i]][[k]][[j]]
				
				#browser()
				
				if (decimal) {
					char.trt.eff <- round(sapply(effCoefList, function(x) x[1]), digits = digits)
					char.trt <- sapply(effCoefList, function(x) 
						ifelse(length(x)>2,  paste(round(x[2:length(x)], digits = digits), collapse = ","), 
										round(x[2], digits = digits)))
				} else {
					char.trt <- sapply(effCoefList, function(x) 
						ifelse(length(x)>2,  paste(attr(fractions(x[2:length(x)]), "fracs"), collapse = ","), 
										attr(fractions(x[2]), "fracs")))
									
					char.trt.eff <- attr(fractions(sapply(effCoefList, function(x) x[1])), "fracs")
				}
		       	
				trt.temp <- c(char.trt, char.trt.eff)
                
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
     trt <- noquote(ifelse(trt == "0", "", trt))
	 
    if (table.legend) {
        Legend <- paste(paste(letters[1:(length(colnames(trt)))], colnames(trt), sep = " = "))
        colnames(trt) <- letters[1:(length(colnames(trt)))]
        trt <- list(trt = trt, Legend = Legend)
    }
    
    
    return(trt)
    
} 
