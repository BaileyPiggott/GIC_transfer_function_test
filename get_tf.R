# function to create transfer functions

get_tf_ind <- function(spec_H, spec_D, spec_Z, spec_A, freq){
  # H, D, Z are three components of magnetic field; A is geomagnetically induced current
  # freq is vector of frequencies at which to calculate TF value
  
  TF_H=vector(mode = 'numeric', length = length(freq)) # initialize empty matrix to store tf
  TF_D=vector(mode = 'numeric', length = length(freq)) # initialize empty matrix to store tf
  TF_Z=vector(mode = 'numeric', length = length(freq)) # initialize empty matrix to store tf
    
  # create vectors from geomag data; each column vector is value of all blocks at one frequency
  # rows correspond to blocks, and columns correspond to one frequency
  F_H <- t(as.matrix(spec_H_block)) # sampling at every frequency
  F_D <- t(as.matrix(spec_D_block))
  F_Z <- t(as.matrix(spec_Z_block))
  F_A <- t(as.matrix(spec_A_block))
  
  
  #calculate and plot transfer function for each component (M)
  for(i in 1:length(freq)){
    
    F_block <- as.matrix(data.frame(F_H[,i], F_D[,i], F_Z[,i]))
    
    M_test[i,] <- ginv(F_block) %*% F_A[,i]
    
    TF_H[i] <- ginv(F_H[,i]) %*% F_A[,i]
    TF_D[i] <- ginv(F_D[,i]) %*% F_A[,i]
    TF_Z[i] <- ginv(F_Z[,i]) %*% F_A[,i]
  }

  # each column is one component of the transfer function at all sampled frequency
  tf <- data.frame('freq' = freq, 'H' = TF_H, 'D' = TF_D, 'Z' = TF_Z) %>% 
    gather(tf, val, 2:4)
  
  
  return(tf)
}