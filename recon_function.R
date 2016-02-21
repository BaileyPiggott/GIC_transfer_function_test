## reconstruction function 
# convert predictions from the frequency domain back to time domain


prediction_fmax <- function(ftest, k, cutoff){# fmax function
  #needed for reconstruction
  
  Fval <- ftest
  
  fMaxInd <- which(Fval > qf(cutoff, 2, 2*k))
  maxes <- c()
  
  if (length(fMaxInd) == 0){
    return(maxes)
  }
  
  for (i in 1:length(fMaxInd)){
    if (fMaxInd[i] == 1 || fMaxInd[i] == length(Fval)){
      next
    }
    
    if (Fval[fMaxInd[i]] > Fval[fMaxInd[i]-1] && 
          Fval[fMaxInd[i]] > Fval[fMaxInd[i]+1]){
      maxes <- c(maxes, fMaxInd[i])
    }
  }
  
  maxes
}

reconstruct <- function(data, block_N, k, nw){ # main reconstruction function
  #need length of time domain block, nw, and k to generate slepians
  
  #complex mean values---------------
  #generate the same slepian tapers used in previous spectral estimates:
  slep <- dpss(n = block_N, k = k, nw = nw, returnEigenvalues = FALSE)$v 
  
  U_kzero <- mvfft(slep)[1, ] # You only want the zeroeth frequency
  
  if (k >= 2){
    U_kzero[seq(2,k,2)] <- 0  # only want the even tapers... see P&W
  }
  ssqU_kzero <- sum(U_kzero^2)
  
  all_pred_cmv <- (data %*% U_kzero) / ssqU_kzero #all cmv's of prediction
  
  #Harmonic F test ----------------   
  pred_ftest <- rep(0, nrow(data)) #initialize vector of f test results
  
  for(i in 1:nrow(data)){ # formula from Thomson, 1982 (equation 13.10)
    denom = 0;
    
    for(j in 1:k){
      denom <- denom + abs(data[i,j] - all_pred_cmv[i]*U_kzero[j])^2
    }
    pred_ftest[i] <- Re(((k-1)*ssqU_kzero*abs(all_pred_cmv[i])^2)/denom)
  }
  
  pred_freqIndex <- prediction_fmax(pred_ftest, k, 0.9) #get indices of significant frequencies
  
  #reconstruct prediction----------
  
  pred_cmv <- rep(0, nrow(data)) #initialize vector of significant cmv's
  pred_cmv[pred_freqIndex] <- all_pred_cmv[pred_freqIndex] #only take frequencies from f test
  pred_cmv[length(all_pred_cmv)-pred_freqIndex+2] <- Conj(all_pred_cmv[pred_freqIndex])
  
  pred_recon <- data.frame(Re(fft(pred_cmv, inverse=TRUE)[1:block_N])) #inverse fourier transform
  
  colnames(pred_recon)<- 'a'
  
  return(pred_recon)
  
}
