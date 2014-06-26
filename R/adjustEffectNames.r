#' Adjust the Effects Names
#'
#' Adjust for appropriate syntax describing the effects matching the 
#' structual formula.
#' 
#' @param effectNames a vector of charactor containing the labels of the
#' treatment or block terms in the model generated by the \code{\link{terms}}.
#' @param effectsMatrix a matrix of variables by terms showing which variables
#' appear in which terms generated by the \code{\link{terms}}.
#' @return a vector of charactor containing the labels of the terms
#' in the model with appropriate syntax describing the effects. 
#' @author Kevin
#' @seealso \code{\link{terms}}
#' @examples
#' 
#' str = "A*(B/E/C)*D"
#' effectsMatrix= attr(terms(as.formula(paste("~", str)), keep.order = TRUE) , "factor")
#' effectNames =  attr(terms(as.formula(paste("~", str)), keep.order = TRUE) , "term.labels")
#' 
#' adjustEffectNames(effectsMatrix, effectNames) 
#' 

adjustEffectNames = function(effectsMatrix, effectNames){
  
   
  nEffect = length(effectNames)
  effectsMatrix = rbind(effectsMatrix,0)
  
  for(j in 1:ncol(effectsMatrix)){
    newEffectNames = ""
    for(i in 1:(nrow(effectsMatrix)-1)){
      
      if(effectsMatrix[i,j] == 0 && effectsMatrix[i+1,j] == 1 &&  (i+1) == nrow(effectsMatrix)){
        newEffectNames = paste(newEffectNames, rownames(effectsMatrix)[i+1], sep = "")  
      } else if(effectsMatrix[i,j] == 0){
        next
      } else if(i != 1 && effectsMatrix[i,j] == 1 &&    effectsMatrix[i-1,j ] == 2 && all(effectsMatrix[(i+1):nrow(effectsMatrix) , j] == 0) ){
        next 
        
      } else if(i != 1 && effectsMatrix[i,j] == 1 &&    effectsMatrix[i-1,j ] == 2 &&   (i+1) == nrow(effectsMatrix) ){
        newEffectNames = paste(newEffectNames, "*", rownames(effectsMatrix)[i+1], sep = "")  
        
        
      } else if(i != 1 && effectsMatrix[i,j] == 1 &&    effectsMatrix[i- 1,j ] == 2 &&   any(effectsMatrix[(i+1):nrow(effectsMatrix), j ] == 1)  ){
        newEffectNames = paste(newEffectNames, "*", sep = "")  
        
      } else if(effectsMatrix[i,j] == 1 && effectsMatrix[i+1,j] == 1 &&   i==1 && (i+1) == nrow(effectsMatrix)){
        newEffectNames = paste(newEffectNames, rownames(effectsMatrix)[i], "*", rownames(effectsMatrix)[i+1], sep = "")
      }else if(effectsMatrix[i,j] == 1 && effectsMatrix[i+1,j] == 1 &&   (i+1) == nrow(effectsMatrix)){
        newEffectNames = paste(newEffectNames, rownames(effectsMatrix)[i+1], sep = "")
      } else if(effectsMatrix[i,j] == 1 &&   any(effectsMatrix[(i+1):nrow(effectsMatrix), j ] == 1)  ){
        newEffectNames = paste(newEffectNames, rownames(effectsMatrix)[i], "*", sep = "")  
        
      } else if(effectsMatrix[i,j] == 1 &&   (all(effectsMatrix[(i+1):nrow(effectsMatrix) , j] == 0) )){
        newEffectNames = paste(newEffectNames, rownames(effectsMatrix)[i], sep = "")  
        
      } else if(effectsMatrix[i,j] == 1 &&   any(effectsMatrix[(i+1):nrow(effectsMatrix), j ] == 2)  ){
        newEffectNames = paste(newEffectNames, rownames(effectsMatrix)[i], "*", "(", sep = "")  
        
      }else if(effectsMatrix[i,j] == 1 &&    effectsMatrix[i- 1,j ] == 1 ){
        next 
        
      }  else if (effectsMatrix[i,j] == 2 &&   effectsMatrix[i+ 1,j ] == 1){
        newEffectNames = paste(newEffectNames, rownames(effectsMatrix)[i], "(",  rownames(effectsMatrix)[i + 1],")", sep = "")  
        
      } else if (effectsMatrix[i,j] == 2 &&   effectsMatrix[i+ 1,j ] == 2){
        newEffectNames = paste(newEffectNames, rownames(effectsMatrix)[i], sep = "")   
        
      }
      
      #print(newEffectNames)
    } 
    
    effectNames[j] = newEffectNames
  } 
  
  return(effectNames)
  
}
