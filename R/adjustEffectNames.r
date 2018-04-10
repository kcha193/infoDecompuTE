adjustEffectNames = function(effectsMatrix, effectNames) {
  nEffect = length(effectNames)
  
  for (j in 1:ncol(effectsMatrix)) {
    newEffectNames = ""
    
    temp = effectsMatrix[, j]
    
    temp = c(temp[temp != 0], 0)
   
    for (i in 1:(length(temp))) {
      if (temp[i] == 1 &&   temp[i + 1] == 1) {
        newEffectNames = paste(newEffectNames, names(temp)[i], "*", sep = "")
        
      } else if (temp[i] == 1 &&   temp[i + 1] == 2) {
        newEffectNames = paste(newEffectNames, names(temp)[i], "*", sep = "")
        
      }  else if (i != 1 && temp[i] == 1 && temp[i - 1] == 2) {
        newEffectNames = paste(newEffectNames, names(temp)[i], "]", sep = "")
        
      } else if (temp[i] == 1 && temp[i + 1] == 0) {
        newEffectNames = paste(newEffectNames, names(temp)[i], sep = "")
      } else  if (temp[i] == 2 && temp[(i + 1):length(temp)][1] == 1) {
        newEffectNames = paste(newEffectNames, names(temp)[i], "[", sep = "")
        
      } else  if (temp[i] == 2 && temp[i + 1] == 2) {
        newEffectNames = paste(unique(unlist(strsplit(c(newEffectNames, names(temp)[i]), "[[:punct:]]"))), ".", sep = "")
        
      }
      
      #print(newEffectNames)
    }
    
    effectNames[j] = newEffectNames
  }
  
  return(effectNames)
  
}



# first to break breakets








