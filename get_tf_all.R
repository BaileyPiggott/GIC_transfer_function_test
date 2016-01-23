# function to create transfer function dependent on all three variables

get_tf_all <- function(spec_H, spec_D, spec_Z, spec_A, freq){

  # H, D, Z are three components of magnetic field; A is geomagnetically induced current
  # freq is vector of frequencies at which to calculate TF value
  tf <- matrix(ncol = 3 , nrow = length(freq)) # initialize empty matrix for transfer function
  
  # create vectors from geomag data; each column vector is value of all blocks at one frequency
  # rows correspond to blocks, and columns correspond to one frequency  
  F_H <- matrix(ncol = length(freq), nrow = (ncol(spec_H))) # initialize empty matrix for sampling at given frequencies
  F_D <- matrix(ncol = length(freq), nrow = (ncol(spec_D))) # initialize empty matrix for sampling at given frequencies
  F_Z <- matrix(ncol = length(freq), nrow = (ncol(spec_Z))) # initialize empty matrix for sampling at given frequencies
  F_A <- matrix(ncol = length(freq), nrow = (ncol(spec_A))) # initialize empty matrix for sampling at given frequencies
  
  
  for(j in 1:length(freq)){
    
    F_H[,j] <- t(spec_H[freq[j],]) #take row corresponding to jth chosen frequency
    F_D[,j] <- t(spec_D[freq[j],]) #take row corresponding to jth chosen frequency
    F_Z[,j] <- t(spec_Z[freq[j],]) #take row corresponding to jth chosen frequency
    F_A[,j] <- t(spec_A[freq[j],]) #take row corresponding to jth chosen frequency
    
  }
 
  #calculate transfer function for each component
  
  for(i in 1:length(freq)){
    
    F_block <- as.matrix(data.frame(F_H[,i], F_D[,i], F_Z[,i]))
    
    tf[i,] <- ginv(F_block) %*% F_A[,i] # each column is the three component tf at given freq

  }
  
  # create data frames of tf magnitude and phase ----------------

  
  return(tf)
}

