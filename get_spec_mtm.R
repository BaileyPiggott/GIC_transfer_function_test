# function to block data and then generate spectral estimate 

get_spec_mtm <- function(data, nw, k, block_length){
  #input data is a vector; N is block length
  #nw and k are the timebandwidth parameter and the number of tapers, respectively
  
  N <- block_length
  blocks <- matrix(ncol = length(data)/N, nrow = N) #initialize empty data blocks
  
  for(j in 1:ncol(blocks)){
    blocks[,j]<- data[(1+(j-1)*N):(j*N)] #each column is a block of data of length N
  }
  
  # get spectral estimate for each column(block) -------
  
  nFFT <- 2^(floor(log2(N)) + 3) #for zero-padding 
  spec <- matrix(nrow = nFFT, ncol = k * ncol(blocks)) #initialize empty matrix for spectral estimate
  
  for(i in 1:(ncol(blocks))){
    
    data <- blocks[,i]
    
    # generate the slepians (the $v will just return the matrix)
    slep <- dpss(n=N, k=k, nw=nw, returnEigenvalues=FALSE)$v
    
    # dummy container to hold the tapered data + extra zeros
    y_k.tmp <- matrix(0, nrow=nFFT, ncol=k)
    
    # multiplies each column of your slepian matrix by the data
    y_k.tmp[1:N, ] <- apply(slep, 2, "*", data)
    
    # This will be a matrix with k columns and nFFT rows - each column is one eigencoefficient
    y_k <- mvfft(y_k.tmp)
    
    spec[,((i-1)*k+1):(i*k)] <- y_k
    
  }
  return(spec)  
}
