


des1 <- local({
  
  Blocks <- factor(rep(1:6, each=3))
  Plots <- factor(rep(1:3, 6))
  A <- factor(rep(1:3, each = 3, times = 2))
  B <- factor(rep(1:3, 6))
  

  data.frame(Blocks, Plots, A, B)
  
}) 
des1


design.df=des1
blk.str = "Blocks/Plots"
trt.str="A*B"

rT <- stats::terms(stats::as.formula(paste("~", blk.str, sep = "")), keep.order = TRUE)  #random terms

rT.terms <- attr(rT, "term.labels")

fT <- stats::terms(stats::as.formula(paste("~", trt.str, sep = "")), keep.order = TRUE)  #fixed terms

#Block structures
Z <- makeBlkDesMat(design.df, rev(rT.terms))

Q <- makeOrthProjectors(Z)
Q

#TReatment structures
trtTerm <- attr(fT, "term.labels")
effectsMatrix <- attr(fT, "factor")

C <- makeContrMat(design.df = design.df, effectNames = trtTerm, 
                  effectsMatrix = effectsMatrix, contr.vec = NA)
C

X <- makeOverDesMat(design.df = design.df, effectNames = trtTerm)
X

Replist <- getTrtRep(design.df, trtTerm)
Rep <- Replist$Rep

C$A


##################################################################
fractions(Q$Blocks)

fractions(Q$`Blocks:Plots`)


#############################################################################
#Between Blocks 
fractions(A_block <- t(X) %*% Q$Blocks %*% X)

fractions(eigen(C$A %*% A_block)$va)
fractions(eigen(C$B %*% A_block)$va)
fractions(eigen(C$`A:B` %*% A_block)$va)


#Between Plots Within Blocks 
fractions(A_Plots <- t(X) %*% Q$`Blocks:Plots` %*% X)

fractions(eigen(C$A %*% A_Plots )$va)
fractions(eigen(C$B %*% A_Plots )$va)
fractions(eigen(C$`A:B` %*% A_Plots )$va)

#############################################################################
(rep_A = diag(2, 9,9))
(rep_B = diag(2, 9,9))
(rep_AB = diag(2, 9,9))




################################################
# inversed and squares of the matrix for the treatment replication for each treatment effect
fractions(rep_A_inv <- diag(1/sqrt(2), 9,9))
fractions(rep_B_inv <- diag(1/sqrt(2), 9,9))
fractions(rep_AB_inv <- diag(1/sqrt(2), 9,9))

#############################################################################
#Between Blocks 
A_block <- t(X) %*% Q$Blocks %*% X

fractions(eigen(rep_A_inv %*% C$A %*% A_block %*% rep_A_inv)$va)
fractions(eigen(rep_B %*% C$B %*% A_block %*% rep_B)$va)
fractions(eigen(rep_AB %*% C$`A:B` %*% A_block %*% rep_AB)$va)


#Between Plots Within Blocks 
A_Plots <- t(X) %*% Q$`Blocks:Plots` %*% X

fractions(eigen(rep_A %*%C$A %*% A_Plots %*% rep_A)$va)
fractions(eigen(rep_B %*%C$B %*% A_Plots %*% rep_B)$va)
fractions(eigen(rep_AB %*% C$`A:B` %*% A_Plots %*% rep_AB)$va)

####################################################################






des1 <- local({
  
  Run <- factor(rep(1:6, each=4))
  Animal <- factor(c(1,6,5,4,3,2,
                     2,1,6,5,4,3,
                     3,2,1,6,5,4,
                     4,3,2,1,6,5))
  Tag <- factor(rep(114:117, times=6))
  Trt <- factor(rep(LETTERS[1:2], times=3)[Animal])
  Time <- factor(rep(c(1:3,3), each=6))
  Sample <- factor(rep(c(1:3,3), each=6))
  Subsample <- factor(rep(c(1,1,1,2), times=6))
  
  data.frame(Run, Tag, Animal, Trt, Sample, Time, Subsample)
  
}) 
des1


library(infoDecompuTE)

# ANOVA for Phase 1 design
des1Ph1 <- with(des1, des1[,c("Animal", "Sample", "Subsample",
                              "Tag", "Trt", "Time")])

design.df <- with(des1Ph1, des1Ph1[order(Animal, Time, Subsample),])


summaryAovOnePhase(design.df=design.df,
                   blk.str = "Animal/Sample", 
                   trt.str="Trt*Time")

design.df=design.df
blk.str = "Animal/Sample"
trt.str="Trt*Time"

rT <- stats::terms(stats::as.formula(paste("~", blk.str, sep = "")), keep.order = TRUE)  #random terms

rT.terms <- attr(rT, "term.labels")

fT <- stats::terms(stats::as.formula(paste("~", trt.str, sep = "")), keep.order = TRUE)  #fixed terms

#Block structures
Z <- makeBlkDesMat(design.df, rev(rT.terms))

Q <- makeOrthProjectors(Z)
Q

#TReatment structures
trtTerm <- attr(fT, "term.labels")
effectsMatrix <- attr(fT, "factor")

C <- makeContrMat(design.df = design.df, effectNames = trtTerm, 
                  effectsMatrix = effectsMatrix, contr.vec = NA)
C

X <- makeOverDesMat(design.df = design.df, effectNames = trtTerm)
X

#######################################################
# Nelder 

eigen(C$Trt %*% t(X) %*% Q$Animal %*% X  %*% C$Trt )

eigen(C$Time %*% t(X) %*% Q$Animal %*% X  %*% C$Time  )

eigen(C$`Trt:Time` %*% t(X) %*% Q$Animal %*% X  %*% C$`Trt:Time`  )


###########################
# John and Williams

A <- t(X) %*% Q$Animal %*% X

fractions(C$Trt)
fractions(C$Trt %*% ginv(A) %*% A)


C$Time
C$Time %*% ginv(A) %*% A

fractions(C$`Trt:Time`)
fractions(C$`Trt:Time` %*% ginv(A) %*% A)

C$Trt %*% A 
A %*% C$Trt  

C$`Trt:Time` %*% A 
A %*% C$`Trt:Time`  

######################################################

A

identityMat(6) - K(6)


tr(C$Trt %*% t(C$Trt))/
  tr(C$`Trt` %*% ginv(A) %*% t(C$`Trt`) )

C$`Trt:Time` %*% ginv(A) %*% A

tr(C$`Trt:Time` %*% t(C$`Trt:Time`))/
  tr(diag(sqrt(Rep[,3])) %*% C$`Trt:Time` %*% ginv(A) %*% t(C$`Trt:Time`) %*%  diag(sqrt(Rep[,3])))
  


Astar <- r.adjust$`Trt:Time` %*% t(X) %*% Q$Animal %*% X %*% r.adjust$`Trt:Time` 
eigen(Astar)


A <- t(X) %*% Q$Animal %*% X
Omaga <- ginv(A)

fractions(C$Trt %*% Omaga %*% C$`Trt:Time`)

fractions(C$Trt %*% Omaga %*% C$`Time`)

fractions(C$`Trt:Time` %*% Omaga %*% C$`Time`)

A <- t(X) %*% Q$`Animal:Sample` %*% X
Omaga <- ginv(A)

fractions(C$Trt %*% Omaga %*% C$`Trt:Time`)

fractions(C$Trt %*% Omaga %*% C$`Time`)

fractions(C$`Trt:Time` %*% Omaga %*% C$`Time`)

eigen(C$Trt %*% A)

eigen(C$`Trt:Time` %*% A)


eigen(A %*% C$`Trt:Time`)



eigen( r.adjust$`Trt:Time` %*% A %*% C$`Trt:Time` %*% r.adjust$`Trt:Time` )


eigen( r.adjust$`Trt:Time` %*%  C$`Trt:Time` %*%A %*%  C$`Trt:Time` %*% r.adjust$`Trt:Time` )




C$`Trt:Time` %*% A 

A %*% C$`Trt:Time`

C$`Trt:Time` %*% A %*% C$`Trt:Time`


eigen( r.adjust$`Trt:Time` %*% C$`Trt:Time` %*% A %*% r.adjust$`Trt:Time` )



eigen( r.adjust$`Trt:Time` %*% C$`Trt:Time` %*% A %*% r.adjust$`Trt:Time` )




eigen(C$`Trt:Time` %*% t(X) %*% Q$`Animal:Sample` %*% X  %*% C$`Trt:Time`  )


Replist <- getTrtRep(design.df, trtTerm)
Rep <- Replist$Rep

r.adjust <- vector(mode = "list", 3)
names(r.adjust) <- names(C)

r.adjust$Trt <- solve(sqrt(diag(Rep[, 1])))
r.adjust$Time <- solve(sqrt(diag(Rep[, 2])))
r.adjust$`Trt:Time` <- solve(sqrt(diag(Rep[, 3])))



################################################################
# Prject into Trt vector subspace in the BEtween Animals stratum 

A_trt =  t(X) %*% Q$Animal %*% X 


eigen(t(X) %*% Q$Animal %*% X )

fractions( t(X) %*% Q$Animal %*% X )

r.adjust$`Trt:Time` %*% t(X) %*% Q$e %*% X %*% r.adjust$`Trt:Time`



eigen(C$Inter %*% t(X) %*% Q$Animal %*% X %*% C$Inter)
eigen(C$Trt %*% t(X) %*% Q$Animal %*% X %*% C$Trt)


fractions(A_trt)

A_trt_ginv <- ginv(A_trt)

fractions(A_trt_ginv)

Q_A_trt <- 
  Q$Animal %*% X %*% ginv(C$Trt  %*% A_trt_ginv %*% C$Trt) %*% t(X) %*% Q$Animal


fractions(Q_A_trt)


eigen(r.adjust$Trt %*%  C$Trt %*% A_trt %*% r.adjust$Trt)

eigen(C$Trt %*% A_trt  %*%  C$Time %*% A_trt)


fractions(C$Trt %*%  t(X) %*% X %*%  C$Trt )

fractions(C$Trt %*%  t(X) %*% X %*%  C$Time)

fractions(C$Trt %*%  t(X) %*% X %*%  C$`Trt:Time` )

fractions(C$Trt  %*% A_trt_ginv %*% C$Trt)

fractions(C$Trt  %*% A_trt_ginv %*% C$Trt*4)

fractions(C$Trt  %*% A_trt_ginv %*% C$Time)
fractions(C$Trt  %*% A_trt_ginv %*% C$`Trt:Time`)









eigen(C$Trt %*% A_trt  %*%  C$`Trt:Time` %*% A_trt)



fractions(C$Trt %*% t(X) %*% Q_A_trt %*% X %*% C$Trt)


fractions(r.adjust$Trt %*% C$Trt %*% t(X) %*% Q_A_trt %*% X %*% C$Trt %*% r.adjust$Trt)

eigen(r.adjust$Trt %*% C$Trt %*% t(X) %*% Q_A_trt %*% X %*% C$Trt %*% r.adjust$Trt)$va
eigen(r.adjust$Trt %*% C$Trt %*% t(X) %*% Q_A_trt %*% X %*% C$Trt %*% r.adjust$Trt)$vec



fractions(C$Time %*% t(X) %*% Q_A_trt %*% X %*% C$Time)


fractions(C$`Trt:Time` %*% t(X) %*% Q_A_trt %*% X %*% C$`Trt:Time`)



C$`Time` %*% t(X) %*% (Q$Animal - Q_A_trt) %*% X %*% C$Time


C$`Trt:Time` %*% t(X) %*% (Q$Animal - Q_A_trt) %*% X %*% C$`Trt:Time`







#Chris's design 1

Blocks <-  factor(rep(1:6, rep(3,6)))
Plots <- factor(rep(1:3,6))
A <-  factor(rep(c(1,2,1, 2,3,2, 3,1,3), 2))
B <-  factor(rep(c(1,2), rep(9,2)))
Test4.df <- data.frame(Blocks, Plots, A, B)


summaryAovOnePhase(design.df=Test4.df,
                   blk.str = "Blocks/Plots", 
                   trt.str="A*B")

# $Fixed$EF
#                     eff.A eff.B eff.A$B
# Between Blocks                           
#   A                  1/3                
#   B                        1            
#   A*B                            1/3    
# Between Blocks(Plots)                    
#   A                  2/3                
#   A*B                            2/3    


design.df=Test4.df

table(Test4.df$A, Test4.df$B)


blk.str = "Blocks/Plots"
trt.str="A*B"

rT <- stats::terms(stats::as.formula(paste("~", blk.str, sep = "")), keep.order = TRUE)  #random terms

rT.terms <- attr(rT, "term.labels")

fT <- stats::terms(stats::as.formula(paste("~", trt.str, sep = "")), keep.order = TRUE)  #fixed terms

#Block structures
Z <- makeBlkDesMat(design.df, rev(rT.terms))

Q <- makeOrthProjectors(Z)
Q

names(Q)

#Treatment structures
trtTerm <- attr(fT, "term.labels")
effectsMatrix <- attr(fT, "factor")

C <- makeContrMat(design.df = design.df, effectNames = trtTerm, 
                  effectsMatrix = effectsMatrix, contr.vec = NA)
C

X <- makeOverDesMat(design.df = design.df, effectNames = trtTerm)
X

Replist <- getTrtRep(design.df, trtTerm)
Rep <- Replist$Rep

r.adjust <- vector(mode = "list", 3)
names(r.adjust) <- names(C)

r.adjust$A <- solve(sqrt(diag(Rep[, 1])))
r.adjust$B <- solve(sqrt(diag(Rep[, 2])))
r.adjust$`A*B` <- solve(sqrt(diag(Rep[, 3])))


#########################################################################
#Project Y onto Treatment A in the Block stratum 

A_trt <-  C$A %*% t(X) %*% Q$Blocks %*% X %*% C$A

A_trt_ginv <- ginv(A_trt)

Q_Block_A <- 
  Q$Blocks %*% X %*% C$A  %*% A_trt_ginv %*% C$A %*% t(X) %*% Q$Blocks


fractions(Q_Block_A)


fractions(C$A %*% t(X) %*% Q_Block_A %*% X %*% C$A)

eigen(r.adjust$A %*% C$A %*% t(X) %*% Q_Block_A %*% X %*% C$A %*% r.adjust$A)$va
eigen(r.adjust$B %*% C$A %*% t(X) %*% Q_Block_A %*% X %*% C$A %*% r.adjust$A)$va
eigen(r.adjust$`A*B` %*% C$A %*% t(X) %*% Q_Block_A %*% X %*% C$`A*B` %*% r.adjust$`A*B`)$va

#########################################################################
#Sweep the  Q_Block_A from Q$Blocks

Q_Block_noA <- Q$Blocks - Q_Block_A

A_trt <-  C$B %*% t(X) %*% Q_Block_noA %*% X %*% C$B

A_trt_ginv <- ginv(A_trt)


#########################################################################
#Project Y onto Treatment B in the Block stratum 

Q_Block_B <- 
  Q$Blocks %*% X %*% C$B  %*% A_trt_ginv %*% C$B %*% t(X) %*% Q$Blocks


fractions(Q_Block_B)

eigen(r.adjust$A %*% C$A %*% t(X) %*% Q_Block_B %*% X %*% C$A %*% r.adjust$A)$va
eigen(r.adjust$B %*% C$B %*% t(X) %*% Q_Block_B %*% X %*% C$B %*% r.adjust$B)$va
eigen(r.adjust$`A*B` %*% C$`A:B` %*% t(X) %*% Q_Block_B %*% X %*% C$`A:B` %*% r.adjust$`A*B`)$va



#########################################################################
#Sweep the  Q_Block_B from Q_Block_noA

Q_Block_noB <- Q_Block_noA - Q_Block_B

A_trt <-  C$`A:B` %*% t(X) %*% Q_Block_noB %*% X %*% C$`A:B`

A_trt_ginv <- ginv(A_trt)


#########################################################################
#Project Y onto Treatment A*B in the Block stratum 

Q_Block_AB <- 
  Q$Blocks %*% X %*% C$`A:B`  %*% A_trt_ginv %*% C$`A:B` %*% t(X) %*% Q$Blocks


fractions(Q_Block_B)

fractions(eigen(r.adjust$A %*% C$A %*% t(X) %*% Q_Block_AB %*% X %*% C$A %*% r.adjust$A)$va)
fractions(eigen(r.adjust$B %*% C$B %*% t(X) %*% Q_Block_AB %*% X %*% C$B %*% r.adjust$B)$va)
fractions(eigen(r.adjust$`A*B` %*% C$`A:B` %*% t(X) %*% Q_Block_AB %*% X %*% C$`A:B` %*% r.adjust$`A*B`)$va)

#########################################################################
#Project Y onto Treatment A in the Plots Within Block stratum 

A_trt <-  C$A %*% t(X) %*% Q$`Blocks:Plots` %*% X %*% C$A

A_trt_ginv <- ginv(A_trt)

Q_Plots_A <- 
  Q$`Blocks:Plots` %*% X %*% C$A  %*% A_trt_ginv %*% C$A %*% t(X) %*% Q$`Blocks:Plots`


fractions(eigen(r.adjust$A %*% C$A %*% t(X) %*% Q_Plots_A %*% X %*% C$A %*% r.adjust$A)$va)
fractions(eigen(r.adjust$B %*% C$B %*% t(X) %*% Q_Plots_A %*% X %*% C$B %*% r.adjust$B)$va)
fractions(eigen(r.adjust$`A*B` %*% C$`A:B` %*% t(X) %*% Q_Plots_A %*% X %*% C$`A:B` %*% r.adjust$`A*B`)$va)


#########################################################################
#Sweep the  Q_Block_A from Q$Blocks

Q_Plots_noA <- Q$`Blocks:Plots` - Q_Plots_A

A_trt <-  C$B %*% t(X) %*% Q_Block_noA %*% X %*% C$B

A_trt_ginv <- ginv(A_trt)


#########################################################################
#Project Y onto Treatment B in the Plots Within Block stratum 

Q_Plots_B <- 
  Q$`Blocks:Plots` %*% X %*% C$B  %*% A_trt_ginv %*% C$B %*% t(X) %*% Q$`Blocks:Plots`


fractions(Q_Block_B)

eigen(r.adjust$A %*% C$A %*% t(X) %*% Q_Block_B %*% X %*% C$A %*% r.adjust$A)$va
eigen(r.adjust$B %*% C$B %*% t(X) %*% Q_Block_B %*% X %*% C$B %*% r.adjust$B)$va
eigen(r.adjust$`A*B` %*% C$`A:B` %*% t(X) %*% Q_Block_B %*% X %*% C$`A:B` %*% r.adjust$`A*B`)$va



#########################################################################
#Sweep the  Q_Block_B from Q_Block_noA

Q_Plots_noB <- Q_Plots_noA - Q_Plots_B

A_trt <-  C$`A:B` %*% t(X) %*% Q_Plots_noB %*% X %*% C$`A:B`

A_trt_ginv <- ginv(A_trt)


#########################################################################
#Project Y onto Treatment A*B in the Block stratum 

Q_Plots_AB <- 
  Q$`Blocks:Plots` %*% X %*% C$`A:B`  %*% A_trt_ginv %*% C$`A:B` %*% t(X) %*% Q$`Blocks:Plots`


fractions(Q_Plots_AB)

fractions(eigen(r.adjust$A %*% C$A %*% t(X) %*% Q_Plots_AB %*% X %*% C$A %*% r.adjust$A)$va)
fractions(eigen(r.adjust$B %*% C$B %*% t(X) %*% Q_Plots_AB %*% X %*% C$B %*% r.adjust$B)$va)
fractions(eigen(r.adjust$`A*B` %*% C$`A:B` %*% t(X) %*% Q_Plots_AB %*% X %*% C$`A:B` %*% r.adjust$`A*B`)$va)


########################################################################################
#Chris's second design 

chris2 <- read.csv("Chris.csv")


summaryAovOnePhase(design.df=chris2,
                   blk.str = "Blocks/Plots", 
                   trt.str="A")

table(chris2$A, chris2$Plots)

design.df=chris2
blk.str = "Blocks/Plots"
trt.str="A+B"

rT <- stats::terms(stats::as.formula(paste("~", blk.str, sep = "")), keep.order = TRUE)  #random terms

rT.terms <- attr(rT, "term.labels")

#Block structures
Z <- makeBlkDesMat(design.df, rev(rT.terms))

Q <- makeOrthProjectors(Z)
Q

names(Q)

#TReatment structures

fT <- stats::terms(stats::as.formula(paste("~", trt.str, sep = "")), keep.order = TRUE)  #fixed terms

trtTerm <- attr(fT, "term.labels")
effectsMatrix <- attr(fT, "factor")

C <- makeContrMat(design.df = design.df, effectNames = trtTerm, 
                  effectsMatrix = effectsMatrix, contr.vec = NA)
C

X <- makeOverDesMat(design.df = design.df, effectNames = trtTerm)
X

Replist <- getTrtRep(design.df, trtTerm)
Rep <- Replist$Rep

r.adjust <- vector(mode = "list", 3)
names(r.adjust) <- names(C)

r.adjust$A <- solve(sqrt(diag(Rep[, 1])))
r.adjust$B <- solve(sqrt(diag(Rep[, 2])))
r.adjust$`A*B` <- solve(sqrt(diag(Rep[, 3])))


#########################################################################
#Project Y onto Treatment A in the Block stratum 

A_trt <-  C$A %*% t(X) %*% Q$Blocks %*% X %*% C$A

A_trt_ginv <- ginv(A_trt)

Q_Block_A <- 
  Q$Blocks %*% X %*% C$A  %*% A_trt_ginv %*% C$A %*% t(X) %*% Q$Blocks


fractions(Q_Block_A)


fractions(C$A %*% t(X) %*% Q_Block_A %*% X %*% C$A)

eigen(r.adjust$A %*% C$A %*% t(X) %*% Q_Block_A %*% X %*% C$A %*% r.adjust$A)$va
eigen(r.adjust$B %*% C$B %*% t(X) %*% Q_Block_A %*% X %*% C$B %*% r.adjust$B)$va


#########################################################################
#Sweep the  Q_Block_A from Q$Blocks

Q_Block_noA <- Q$Blocks - Q_Block_A

A_trt <-  C$B %*% t(X) %*% Q_Block_noA %*% X %*% C$B

fractions(A_trt) #Stop here!!

A_trt_ginv <- ginv(A_trt)


#########################################################################
#Project Y onto Treatment A in the Plots Within Block stratum 

A_trt <-  C$A %*% t(X) %*% Q$`Blocks:Plots` %*% X %*% C$A

A_trt_ginv <- ginv(A_trt)

Q_Plots_A <- 
  Q$`Blocks:Plots` %*% X %*% C$A  %*% A_trt_ginv %*% C$A %*% t(X) %*% Q$`Blocks:Plots`


fractions(eigen(r.adjust$A %*% C$A %*% t(X) %*% Q_Plots_A %*% X %*% C$A %*% r.adjust$A)$va)
fractions(eigen(r.adjust$B %*% C$B %*% t(X) %*% Q_Plots_A %*% X %*% C$B %*% r.adjust$B)$va)

#########################################################################
#Sweep the  Q_Block_A from Q$Blocks

Q_Plots_noA <- Q$`Blocks:Plots`   - Q_Plots_A

A_trt <-  C$B %*% t(X) %*% Q_Plots_noA %*% X %*% C$B

A_trt_ginv <- ginv(A_trt)


#########################################################################
#Project Y onto Treatment B in the Plots Within Block stratum 

Q_Plots_B <- 
  Q_Plots_noA %*% X %*% C$B  %*% A_trt_ginv %*% C$B %*% t(X) %*% Q_Plots_noA


fractions(Q_Plots_B)

eigen(r.adjust$A %*% C$A %*% t(X) %*% Q_Plots_B %*% X %*% C$A %*% r.adjust$A)$va
eigen(r.adjust$B %*% C$B %*% t(X) %*% Q_Plots_B %*% X %*% C$B %*% r.adjust$B)$va





