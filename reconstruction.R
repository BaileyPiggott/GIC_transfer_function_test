## recontruction of data

#Functions -------------------
# Finds local maxes in the Ftest given a spec.mtm object.
# based on a probability cutoff

#sketchy version of find local Fmax function

test_function <- function(ftest, k, cutoff){
  
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



# Returns the indices.
findLocalFMax <- function(obj, cutoff){
  # Check whether this is a spec.mtm() object, or from my own CMV code.
  if (any(class(obj) == "Ftest")){
    Fval <- obj$mtm$Ftest
    k <- obj$mtm$k
  } else if (any(class(obj) == "driegert.cmv")){
    Fval <- obj$Ftest$Fval
    k <- obj$k
  } else {
    stop("obj needs to be of class 'driegert.cmv' or 'Ftest'.")
  }
  
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

#load data ----------
source('get_tf_all.R')
library(MASS) # need for pseudoinverse

sim_data <- as.data.frame(read.table("MA6.txt"))
h <- sim_data$V2
d <- sim_data$V3
z <- sim_data$V4
a <- sim_data$V5 #current 

x <- 1:length(a)
fit <- lm(a ~ x + I(x^2) + I(x^3))
cubeTrend <- fitted.values(fit)
a2 <- a - cubeTrend # removing a cubic trend (want this for the reconstruction part)

#plot check----------------
## top: plot with cubic fit overlay
## bottom: residuals after removing cubic fit.
par(mfrow=c(2,1), mar=c(4,4,1,1), mgp=c(2.5,1,0))
plot(a, type='l', main = "Original Data with Cubic Trend")
lines(cubeTrend, col='red', lwd=2)
plot(a2, type='l', ylab="Residuals", main = "Residuals") # want 0-mean for reconstruction

#calculate spectrum for actual data ----------------------
require('multitaper')

# calculate the spectrum (returnInternals gives us the complex-mean-values (cmv's))

# cmv's are an average of the eigencoefficients (windowed-fourier-transformed data)
# if you average the squared eigencoefficients, you get the spectrum.

s <- spec.mtm(a2, Ftest=TRUE, returnInternals=TRUE, plot=FALSE)

freqIdx <- findLocalFMax(s, 0.9)# find frequencies of large F-values ("most sinusoidal" frequencies)

# set up the CMV array (spec.mtm() only returns the positive frequencies due to symmetry)
cmv <- rep(0, s$mtm$nFFT)
cmv[freqIdx] <- s$mtm$cmv[freqIdx]
cmv[s$mtm$nFFT-freqIdx+2] <- Conj(s$mtm$cmv[freqIdx])

# the "negative" side of the cmv's actually go in the right side of the array and must be
# conjugated (due to the minus sign on the frequency)

#create transfer function for prediction -----------

sh <- spec.mtm(h, Ftest=TRUE, returnInternals=TRUE, plot=FALSE)
sd <- spec.mtm(d, Ftest=TRUE, returnInternals=TRUE, plot=FALSE)
sz <- spec.mtm(z, Ftest=TRUE, returnInternals=TRUE, plot=FALSE)

freq <- seq(1,nrow(sh$mtm$eigenCoefs),1) # frequencies at which to calculate the transfer function

tf <- get_tf_all(sh$mtm$eigenCoefs, sd$mtm$eigenCoefs, sz$mtm$eigenCoefs, s$mtm$eigenCoefs, freq)

#predict current ------------------

F_test_H <- matrix(ncol = length(freq), nrow = 7) # initialize empty matrix for sampling at given frequencies
F_test_D <- matrix(ncol = length(freq), nrow = 7) # initialize empty matrix for sampling at given frequencies
F_test_Z <- matrix(ncol = length(freq), nrow = 7) # initialize empty matrix for sampling at given frequencies


for(j in 1:length(freq)){  
  F_test_H[,j] <- t(sh$mtm$eigenCoefs[freq[j],]) #take row corresponding to jth chosen frequency
  F_test_D[,j] <- t(sd$mtm$eigenCoefs[freq[j],]) #take row corresponding to jth chosen frequency
  F_test_Z[,j] <- t(sz$mtm$eigenCoefs[freq[j],]) #take row corresponding to jth chosen frequency
}

F_predict_A <- matrix(ncol = length(freq), nrow = 7) # initialize empty matrix for sampling at given frequencies

for(i in 1:length(freq)){
  F_block <- as.matrix(data.frame(F_test_H[,i], F_test_D[,i], F_test_Z[,i]))  
  F_predict_A[,i] <- F_block %*% tf[i,]
}

prediction <- matrix(ncol = 7, nrow = nrow(sh$mtm$eigenCoefs))

for(j in 1:length(freq)){
  prediction[freq[j],] <- t(F_predict_A[,j])
}

#reconstruction of actual data ---------------

recon <- Re(fft(cmv, inverse=TRUE)[1:length(a2)])
# the fourier transform is complex, however since we started with real data, 
# the imaginary part will be zero, so just take the real part.

fullRecon <- recon + cubeTrend

par(mfrow=c(2,1))
plot(a2, type='l', col='grey', lwd=2, main="Reconstruction of Zero Mean Data")
lines(recon, col='blue')
plot(a, type='l', col='grey', lwd=2, main="Reconstruction of Original Data")
lines(fullRecon, col='blue')

#cmv's of prediction -----------

#assume F-test from training data applies to prediction
ftest <- s$mtm$Ftest
k=7
freqIndex_pred <- test_function(ftest, k, 0.9)

#complex mean values:
slep <- dpss(n = length(a), k = k, nw = 4, returnEigenvalues = FALSE)$v #slepian tapers

U_kzero <- mvfft(slep)[1, ] # You only want the zeroeth frequency

if (k >= 2){
  U_kzero[seq(2,k,2)] <- 0  # only want the even tapers... see P&W
}
ssqU_kzero <- sum(U_kzero^2)

pred_cmv <- (prediction %*% U_kzero) / ssqU_kzero

#reconstruct prediction----------

test_cmv <- rep(0, nrow(prediction))
test_cmv[freqIndex_pred] <- pred_cmv[freqIndex_pred] #only take frequencies from f test
test_cmv[length(pred_cmv)-freqIndex_pred+2] <- Conj(pred_cmv[freqIndex_pred])

test_recon <- Re(fft(test_cmv, inverse=TRUE)[1:length(a)])

#plot prediction------------

par(mfrow=c(2,1))
plot(a, type='l', col='grey', lwd=2, main="Reconstruction of Prediction")
lines(test_recon, col='blue')
plot(a, type='l', col='grey', lwd=2, main="add cubic trend")
lines(recon + cubeTrend, col='blue')

