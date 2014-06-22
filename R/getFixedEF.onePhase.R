getFixedEF.onePhase <- function(effFactors, trt.Coef, T, Rep, table.legend, decimal, 
    digits) {
    
    trt <- numeric(length(trt.Coef) + ncol(Rep))
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
            
            char.trt <- attr(fractions(trt.Coef * unlist(effFactors[[i]][[j]])), "fracs")
            char.trt.eff <- attr(fractions(unlist(effFactors[[i]][[j]])), "fracs")
            
            if (decimal) {
                char.trt <- round(trt.Coef * unlist(effFactors[[i]][[j]]), digits = digits)
                char.trt.eff <- round(unlist(effFactors[[i]][[j]]), digits = digits)
            }
            
            trt.temp <- c(char.trt, char.trt.eff)
            
            trt <- rbind(trt, trt.temp)
            rownames(trt)[nrow(trt)] <- paste("  ", names(effFactors[[i]][j]), sep = " ")
        }
        
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
