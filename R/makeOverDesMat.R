makeOverDesMat <- function(design.df, effectNames) {
    
    if (length(effectNames) == 1 && ! any(grepl(":", effectNames))) {
        incident = design.df[, effectNames]
        nLevels = levels(design.df[, effectNames])
    } else if (any(grepl(":", effectNames))) {
        uniqueTrtCols = unique(unlist(strsplit(effectNames, "\\:")))
        
        incident = as.factor(apply(design.df[, uniqueTrtCols], 1, function(x) paste(x, collapse = ".")))
        # nLevels = sort(unique(incident))
        nLevels = sort(levels(interaction(design.df[, uniqueTrtCols])))
        
    } else {
        incident = as.factor(apply(design.df[, effectNames], 1, function(x) paste(x, collapse = ".")))
        nLevels = sort(levels(interaction(design.df[, effectNames])))
    }
    
    N <- matrix(0, nrow = nrow(design.df), ncol = length(nLevels))
    N[cbind(1:nrow(design.df), match(incident, nLevels))] <- 1
    
    return(N)
} 
