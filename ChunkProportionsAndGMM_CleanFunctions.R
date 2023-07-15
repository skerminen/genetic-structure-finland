##########################################
#                                        #
#     CHUNK PROPORTIONS AND GMM          #
#                                        #
#     Sini Kerminen and Matti Pirinen    #
#     22.4.2016                          #
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

# Here we describe the code that is used to asses the probability of an 
# individual belonging to either of two compared populations in
# Kerminen et al. (YYYY).


### Input
# The code takes tree files as an input that are described below. 
# coanc = Coancestry matrix from ChromoPainter program (Lawson et al. 2012)
# refindividuals1 = List of first population's reference individual's ids as 
#                   in coancestry matrix
# refindividuals2 = List of second population's reference individual's ids as 
#                   in coancestry matrix


### Output
# The code generates output in data.frame described below.
# results = Data.frame with columns "ID", "FIXED", "PROPORTION", "LOG", "PROB"
#           ID         = Ids of the individuals from the coancestry matrix. These have
#                        to match with the ids in the refindividuals 1 and 2.
#           FIXED      = Addresses individual's status as reference individual (1 or 2) 
#                        or as test individual (0).
#           PROPORTION = Individual's chunk proportions (See Kerminen et al. (YYYY) )
#           LOG        = Chunk proportions in log scale
#           PROB       = Individual's propability to belong to a population based on
#                        a supervised Gaussian mixture model on the log scale chunk proportion values
#                        using an EM-algorithm.



###### FUNCTIONS #########################

chunkProportions <- function(coancrow, ref1idlist, ref2idlist){
  meanRef1 <- sum( coancrow[ names(coancrow) %in% ref1idlist] ) / length(ref1idlist)  
  meanRef2 <- sum( coancrow[ names(coancrow) %in% ref2idlist ]) / length(ref2idlist)   
  return( meanRef1 / meanRef2 )
}

em.gaussian.mixture<-function(y,ncl,fixed=NULL,tol=1e-6,max.iter=10000){
  nind = length(y)
  mu = seq(-1,1,length.out = ncl)
  sigma = rep(1,ncl)
  w.clu = rep(0.5, ncl)
  w.ind = matrix(1/ncl,nrow=nind,ncol=ncl)
  w.fixed = w.ind
  if(!is.null(fixed)){
    for( i in 1:ncl){w.fixed[fixed == i,] = 0;  w.fixed[fixed == i,i] = 1;}
    w.fixed = w.fixed[fixed > 0,]
  }
  res.loglkhood = c()
  loglkhood = univar.gaussian.loglkhood(y,mu,sigma,w.clu)
  iter=0;converged = FALSE;
  while(!converged){
    iter = iter + 1
    a = w.ind*NA
    for( i in 1:ncl){a[,i] = w.clu[i]*dnorm(y,mu[i],sigma[i],log=FALSE)}
    w.ind = a/rowSums(a) #new weights for individual memberships
    if(!is.null(fixed)) w.ind[fixed > 0,] = w.fixed
    w.clu = apply(w.ind,2,mean) #new weights for clusters
    mu = colSums(w.ind*y)/colSums(w.ind) #new means
    sigma = sqrt(rowSums(t(w.ind)*(t(matrix(y,ncol=ncl,,nrow=nind,byrow=FALSE))-mu)^2)/colSums(w.ind)) #new.sds   
    res.loglkhood = c(res.loglkhood,univar.gaussian.loglkhood(y,mu,sigma,w.clu))
    if(iter > max.iter | (res.loglkhood[iter]-loglkhood) < tol) converged = TRUE
    loglkhood = res.loglkhood[iter]
    print(paste(iter,loglkhood))
  }
  return(list(mu=mu,sigma=sigma,w.clu=w.clu,prob.ind=w.ind,loglkhood=loglkhood))
}

univar.gaussian.loglkhood <- function (y,mu,sigma,w){
  ncl = length(mu)
  a = matrix(NA,ncol=ncl,nrow=length(y))
  for( i in 1:ncl){a[,i] = w[i]*dnorm(y,mu[i],sigma[i],log=FALSE)}
  return(sum(log(rowSums(a))))
}


###### EXAMPLE SCRIPT ####################

# Read input 
coanc <- read.table("examples/Example_data.chunkcounts.out", h = T, row.names = 1) # coancestry matrix
refindividuals1 <- paste("IND", 1:10, sep = "") # Lets use first ten inds from the coanc as pop1 reference individuals
refindividuals2 <- paste("IND", 91:100, sep = "") # and last 10 inds as reference inds from other pop

# Initialize results
results <- as.data.frame( matrix(NA, nrow = nrow(coanc), ncol = 1) )
results[,"ID"] <- colnames(coanc)
results[, "FIXED"] <- 0
results[ results[,"ID"] %in% refindividuals1 , "FIXED"] <- 1
results[ results[,"ID"] %in% refindividuals2 , "FIXED"] <- 2
results <- results[,-1]

# Calculate chunk proportions
results[, "PROPORTION"] <- apply( coanc, 1, chunkProportions, refindividuals1, refindividuals2 )
results[,"LOG"] <- log( results[, "PROPORTION"] )

# GMM
res.em = em.gaussian.mixture( results[,"LOG"], 2, results[,"FIXED"] )
results[, "PROB"] <- res.em$prob.ind[,2]

# Check results
head(results)
table(round(results$PROB, digits = 1))

