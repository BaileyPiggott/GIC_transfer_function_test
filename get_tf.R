# function to create transfer functions dependant only on one given variable

get_tf_ind <- function(spec_H, spec_D, spec_Z, spec_A, freq){
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
  
  
  #calculate and plot transfer function for each component (M)
  for(i in 1:length(freq)){
 
    TF_H[i] <- ginv(F_H[,i]) %*% F_A[,i]
    TF_D[i] <- ginv(F_D[,i]) %*% F_A[,i]
    TF_Z[i] <- ginv(F_Z[,i]) %*% F_A[,i]
  }

  # each column is one component of the transfer function at all sampled frequency
  tf <- data.frame('freq' = freq, 'H' = TF_H, 'D' = TF_D, 'Z' = TF_Z) %>% 
    gather(tf, val, 2:4)
  
  
  return(tf)
}