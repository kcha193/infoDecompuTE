getFixedEF.onePhase <- function(effFactors, trt.Sca, T, Rep, table.legend, decimal, 
    digits) {
    
    trt <- numeric(length(trt.Sca) + ncol(Rep))
    names(trt) <- c(names(T), paste("eff", names(T), sep = "."))
    
    for (i in 1:length(effFactors)) {
        trt <- rbind(trt, character(length = length(T) * 2))
        if (grepl("Within", names(effFactors[i]))) {
            rownames(trt)[nrow(trt)] <- paste(names(effFactors[i]), sep = " ")
        } else {
            rownames(trt)[nrow(trt)] <- paste("Between", names(effFactors[i]), sep = " ")
        }
        
        for (j in 1:length(effFactors[[i]])) {
            if (is.null(effFactors[[i]][[j]])) 
                next
            
					
			effCoefList = effFactors[[i]][[j]]
				
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
            rownames(trt)[nrow(trt)] <- paste("  ", names(effFactors[[i]][j]), sep = " ")
        }
        
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
