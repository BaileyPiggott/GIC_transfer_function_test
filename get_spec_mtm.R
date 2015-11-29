# function to generate spectral estimate of data

get_spec_mtm <- function(data, nw, k, block_length){
  #input data is a vector; N is block length
  # generally, nw = 4, k = 7
  
  N <- block_length
  
  blocks <- matrix(ncol = length(data)/N, nrow = N)
  
  for(j in 1:ncol(blocks)){
    blocks[,j]<- data[(1+(j-1)*N):(j*N)]
  }
  
  # get spectral estimate for each column -------
  
  nFFT <- 2^(floor(log2(N)) + 3) 
  spec <- matrix(nrow = nFFT, ncol = k * ncol(blocks))
  
  for(i in 1:(ncol(blocks))){
    
    # this is for zero-padding - you can change the 3 to whatever you want, but 2 or 3 is usually fine
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