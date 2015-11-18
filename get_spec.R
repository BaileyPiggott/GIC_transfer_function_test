# function to generate spectral estimate of data

get_spec <- function(data, block_length){
#input data is a vector; N is block length
    
  N <- block_length
  
  block <- data.frame(matrix(ncol = (length(data)/N), nrow = N))  #initialize empty matrix for data blocks
  spec_est <- data.frame(matrix(ncol = (length(data)/N), nrow = N))  #initialize empty matrix for spectral estimation
  
  for(j in 1:(length(data)/N)){
    
    #separate into blocks:  
    block[,j]<- data[(1+(j-1)*N):(j*N)]
   
    #spectral estimate of each block:
    spec_est[,j] <- (1/N)*abs(fft(block[,j]))^2
  }

  spec_est <- spec_est[1:(N/2+1),] #spectra are symmetric; only need half
  
  
  return(spec_est)
}