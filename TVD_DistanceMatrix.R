##########################################
#                                        #
#         TVD                            #
#                                        #
#         Sini Kerminen                  #
#         30.3.2017                      #
#                                        #
##########################################

###### LICENCE ###########################

# This code is covered by the MIT license.
# 
# 
# The MIT License
# 
# Copyright (c) 2016 University of Helsinki
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

###### DESCRIPTION #######################

### Input
# The code takes tree files as an input that are described below. 
# coanc = Coancestry matrix from ChromoPainter program (Lawson et al. 2012)
# pops  = Population assignment for individuals in coanc. Example:
#         ID   POP1 POP2 POP3
#         IND1 1    1    1
#         IND2 1    2    2
#         IND3 1    2    3


### Output
# Output is a distance matrix that can be illustrated for example as a dendrogram
# with hclust and accompanied functions.


#########################
## FUNCTIONS

## TVD
# Returns TVD value (Leslie et al. 2015). Arguments: 2 copying vectors
calcTVD <- function(copyvector1, copyvector2){
  
  tvd <- 0.5 * sum( abs( copyvector1 - copyvector2) )
  
  return(tvd)
}

## CopyingVector
# Returns a copy vector for a populations as a list. Arguments: population on coancestry
# matrix, population structure and tvd level.
newCopyVector <- function(population, pops, tvdlevel){
  
  print("New copy vector")
  
  copyvector <- vector(length = tvdlevel)
  
  for(i in 1:tvdlevel){
    donorpopinds <- pops[ pops[, paste("POP", tvdlevel, sep = "")] == i, "ID"]
    
    if(i == 1){
      indcopyvectors <- apply(population, 1, indProportions, donorpopinds)
    }
    copyvector[i] <- mean( apply(population, 1, indProportions, donorpopinds)  )
  }
  
  return(copyvector)
}

indProportions <- function(ind, donors){ sum( ind[ as.character(donors) ] ) / sum(ind)  } 


## TVD matrix
# Returns POPxPOP TVD matrix. Arguments: Copying vector matrix
calcTVDmatrix <- function( copyvectormatrix ){
  
  print("New TVD matrix")
  
  tvdmatrix <- matrix(NA, ncol = nrow(copyvectormatrix), nrow = nrow(copyvectormatrix) )
  
  for(i in 1:(nrow(tvdmatrix)-1) ){
    for(j in (i+1):nrow(copyvectormatrix) )
    
      tvdmatrix[i, j] <- calcTVD(copyvectormatrix[i,], copyvectormatrix[j,])
  }
  
  return(tvdmatrix)
}


## Copying vector matrix for TVD tree
# Returns copying vectors for each population as POPxDonorPOP matrix. Arguments: 
# coancestry matrix, population structure, population and tvd level, and pop list
calcCopyVectorMatrixForTree <- function(coanc, pops, poplevel, tvdlevel, list){
  
  copyVectorMatrix <- matrix(0, nrow = length(list), ncol = tvdlevel)
  
  for(i in 1:length(list)){
    
    #define population here
    popinds <- pops[ pops[, paste("POP", poplevel, sep = "")] %in% list[[i]] , "ID" ]
    population <- coanc[ rownames(coanc) %in% popinds , ]  
    
    # Calculate copy vector for population
    copyVectorMatrix[i, ] <- newCopyVector(population, pops, tvdlevel)
  }
  
  return(copyVectorMatrix)
}

## Merge population lists
mergeLists <- function( list, index1, index2 ){
  
  newlist <- list
  newlist[[index1]] <- append(list[[index1]], list[[index2]])
  newlist[[index2]] <- NULL
  
  return(newlist)
}



#########################
## SCRIPT

# Read input data
coanc <- read.table("ChunkProportionsAndGMM/Example_data.chunkcounts.out", header = T, row.names = 1)
pops <- read.table("ChunkProportionsAndGMM/Example_data_populations.txt", header = T)
poplevel <- 4
tvdlevel <- 4

# Initialize output
distancematrix <- matrix(NA, ncol = poplevel, nrow = poplevel)
diag(distancematrix) <- poplevel
mergedpops <- list()
poplist <-  as.list( 1:poplevel )
mergedpops[[1]] <- poplist


# Calculate distance matrix
scale <- 0
for( round in 1:(poplevel-1)){
  
  print(mergedpops[[round]])
  
  # Calculate copy vector and TVD matrix for new population configuration
  copyvectors <- calcCopyVectorMatrixForTree(coanc, pops, poplevel, tvdlevel, mergedpops[[round]])
  tvdmatrix <- calcTVDmatrix( copyvectors )
  
  # Find a pair with minimum TVD value
  minIndeces <- which(tvdmatrix == min(tvdmatrix, na.rm = T), arr.ind = T)
  
  scale <- scale + min( tvdmatrix, na.rm = T )
  
  ## Merge population lists with minimum TVD values
  mergedpops[[round+1]] <- mergeLists( mergedpops[[round]], minIndeces[1,1], 
                                       minIndeces[1,2] )
  
  ## Fill in the distance matrix (symmetric)
  for(k in 1:length(mergedpops[[round]][[minIndeces[1,1] ]]) ){
    for(l in 1:length(mergedpops[[round]][[minIndeces[1,2] ]]) ){
      
      distancematrix[ mergedpops[[round]][[minIndeces[1,1]]][k], mergedpops[[round]][[minIndeces[1,2]]][l]  ] <- scale
      distancematrix[ mergedpops[[round]][[minIndeces[1,2]]][l], mergedpops[[round]][[minIndeces[1,1]]][k]  ] <- scale
      
    }
  }
  
}

# See a symmetric distance matrix as a result
distancematrix