# function to generate spectral estimate of data

get_spec_mtm <- function(data, block_length){
  #input data is a vector; N is block length
  
  N <- block_length
  
  spec_rows <- 2^ceiling(log2(N))+1 # number of rows is a power of 2 larger than the length of the data 
  spec_est <- data.frame(matrix(ncol = (length(data)/N), nrow = spec_rows))  #initialize empty matrix for spectral estimation
  
  for(j in 1:(length(data)/N)){
    
    #separate into blocks:  
    block <- ts(data[(1+(j-1)*N):(j*N)])
    
    #spectral estimate of each block:
    estimate <- spec.mtm(block, plot = FALSE)
    spec_est[,j] <- estimate[['spec']]
  }
  
  #freq <- data.frame('freq' = estimate[['freq']])
  #spec_est <- bind_cols(freq, spec_est) 
  
  return(spec_est)
}