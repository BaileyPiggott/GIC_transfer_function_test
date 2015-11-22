# function to create transfer function dependent on all three variables

get_tf_all <- function(spec_H, spec_D, spec_Z, spec_A, freq){
  # H, D, Z are three components of magnetic field; A is geomagnetically induced current
  # freq is vector of frequencies at which to calculate TF value
  
  TF_H=vector(mode = 'numeric', length = length(freq)) # initialize empty matrix to store tf
  TF_D=vector(mode = 'numeric', length = length(freq)) # initialize empty matrix to store tf
  TF_Z=vector(mode = 'numeric', length = length(freq)) # initialize empty matrix to store tf
  
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
  
  M_test <- matrix(nrow = 3 , ncol = length(freq)) # initialize empty matrix for transfer function
  
  #calculate and plot transfer function for each component (M)
  
  for(i in 1:length(freq)){
    
    F_block <- as.matrix(data.frame(F_H[,i], F_D[,i], F_Z[,i]))
    
    M_test[,i] <- ginv(F_block) %*% F_A[,i] # each column is the three component tf at given freq

  }
  
  return(M_test)
}

