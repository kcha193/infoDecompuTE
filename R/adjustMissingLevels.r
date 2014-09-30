
adjustMissingLevels = function(design.df, str.for){
	
	if(!grepl("(:)", str.for) && !grepl("(/)", str.for)) 
		return(list(design.df = design.df, str.for = str.for))
	#browser()
	fT = terms(as.formula(paste("~", str.for, sep = "")), keep.order = TRUE)

	trtTerm <- attr(fT, "term.labels")
    effectsMatrix <- attr(fT, "factor")
    
    trtTerm = adjustEffectNames(effectsMatrix, trtTerm)
		
	oldTerm = unlist(strsplit(trtTerm, "\\*"))
	
	oldTerm = gsub("\\)", "", unique(oldTerm[grep("\\(", oldTerm)]))
	
	for(i in oldTerm){
	
		oldTermNames = unique(unlist(strsplit(i, "\\(")))
		
		oldTermNames = oldTermNames[oldTermNames!=""]
		
		bwTermName = drop.levels(interaction(design.df[,unlist(strsplit(oldTermNames[1], "\\."))]))
		
		oldTermNames = oldTermNames[2]
		
		newDes = design.df[,oldTermNames]

		newDes = data.frame(bwTermName = as.factor(as.character(bwTermName)), newDes)
		
		if(any(is.na(newDes[,1]))){		
			temp = as.character(newDes[,1])
			temp[is.na(temp)] = "NA1"
			newDes[,1]= temp
		}
		
		if(any(is.na(newDes[,2]))){		
			temp = as.character(newDes[,2])
			temp[is.na(temp)] = "NA1"
			newDes[,2]= temp
		}

		
		newDes = drop.levels(newDes)
		
		colnames(newDes)[2] = termNames = paste("m", oldTermNames, sep = "") 

		termTableList = apply(newDes, 2, table)
		
		#if(!all(sapply(termTableList, function(x) all(outer(x,x,"=="))))) #unbalanced design
		#	next
		
		tempCheck = tapply(newDes[,2], as.factor(newDes[,1]), function(x) x)

		repeated = FALSE
		
		for(i in 1:length(tempCheck)){

			if(any(outer(tempCheck[[i]], unlist(tempCheck[-i]), "=="))){
				repeated = TRUE
				break
			}

		}
		
		if(repeated) 
			next
		
		termListLength = sapply(termTableList,length)

		#The level length of first factor should not be greater than second 
		#if(termListLength[1] >= termListLength[2]) 			
		#	next
				
		check  = names(termListLength)
		inter = interaction(newDes[,check])
		
		if(length(unique(as.character(inter))) != length(levels(inter))){
			
			shortTerm = names(termListLength)[1]
			
			inter = interaction(newDes[,shortTerm])
		
			levelBased = which(table(newDes[,shortTerm]) == max(table(newDes[,shortTerm])))[1]

			repeatedLevels = drop.levels(newDes[,names(termListLength)[2]][inter == levels(inter)[levelBased]])
			
			for(k in  levels(inter)[-levelBased]){
				temp = drop.levels(newDes[,names(termListLength)[2]][inter == k])
				levels(temp) = levels(repeatedLevels)
				newDes[,names(termListLength)[2]][inter == k] = temp
			}
			
			design.df = cbind(design.df, newDes[,names(termListLength)[2]])
			design.df[design.df=="NA2"] = NA			
			names(design.df)[length(design.df)] = names(termListLength)[2]
			str.for = gsub(oldTermNames, names(termListLength)[2], str.for)
		}
			
			
	}	
	
	design.df = drop.levels(design.df)
	
	return(list(design.df = design.df, str.for = str.for))
	
}

