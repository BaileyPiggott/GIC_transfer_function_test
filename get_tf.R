# function to create transfer functions dependant only on one given variable


get_tf_ind <- function(data_in, data_out, freq){
  
  tf <- vector(mode = 'numeric', length = length(freq))
  
  F_in <- matrix(ncol = length(freq), nrow = ncol(data_in))
  F_out <- matrix(ncol = length(freq), nrow = ncol(data_out))
  
  for(j in 1:length(freq)){
    
    F_in[,j] <- t(data_in[freq[j],]) #take row corresponding to jth chosen frequency
    F_out[,j] <- t(data_out[freq[j],]) #take row corresponding to jth chosen frequency
    
  }
  
  for(i in 1:length(freq)){
    tf[i] <- ginv(F_in[,i]) %*% F_out[,i]
  }
  
  # create data frames of tf magnitude and phase ----------------
  p=1/2*length(tf)
  mag <- abs(tf)
  phase <- atan2(Im(tf), Re(tf))
  
  tf <- data.frame(freq = seq(0,0.5,0.5/p), Magnitude = mag[1:(p+1)], Phase = phase[1:(p+1)])
  tf <- tf %>% gather(type, val, Magnitude:Phase)
  
  return(tf)
}