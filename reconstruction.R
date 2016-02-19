## recontruction of data

#Functions -------------------
# Finds local maxes in the Ftest given a spec.mtm object.
# based on a probability cutoff
# Returns the indices.

#version of find local Fmax function for predictions

prediction_fmax <- function(ftest, k, cutoff){
  
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

#dave's function:
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
train_data <- sim_data[1:(nrow(sim_data)/2), ]
test_data <- sim_data[(nrow(sim_data)/2 +1):nrow(sim_data), ]

h_train <- train_data$V2
d_train <- train_data$V3
z_train <- train_data$V4
a_train <- train_data$V5 #current 

h_test <- test_data$V2
d_test <- test_data$V3
z_test <- test_data$V4
a_test <- test_data$V5 #current 


x <- 1:length(a_train)
fit <- lm(a_train ~ x + I(x^2) + I(x^3))
cubeTrend <- fitted.values(fit)
a2 <- a_train - cubeTrend # removing a cubic trend (want this for the reconstruction part)

#plot check----------------
## top: plot with cubic fit overlay
## bottom: residuals after removing cubic fit.
par(mfrow=c(2,1), mar=c(4,4,1,1), mgp=c(2.5,1,0))
plot(a_train, type='l', main = "Original Data with Cubic Trend")
lines(cubeTrend, col='red', lwd=2)
plot(a2, type='l', ylab="Residuals", main = "Residuals") # want 0-mean for reconstruction

#calculate spectrum for actual data ----------------------
require('multitaper')

# calculate the spectrum (returnInternals gives us the complex-mean-values (cmv's))

s <- spec.mtm(a2, Ftest=TRUE, returnInternals=TRUE, plot=FALSE)

freqIdx <- findLocalFMax(s, 0.9)# find frequencies of large F-values ("most sinusoidal" frequencies)

# set up the CMV array (spec.mtm() only returns the positive frequencies due to symmetry)
cmv <- rep(0, s$mtm$nFFT)
cmv[freqIdx] <- s$mtm$cmv[freqIdx]
cmv[s$mtm$nFFT-freqIdx+2] <- Conj(s$mtm$cmv[freqIdx])
# the "negative" side of the cmv's actually go in the right side of the array and must be
# conjugated (due to the minus sign on the frequency)

#create transfer function for prediction -----------

sh_train <- spec.mtm(h_train, Ftest=TRUE, returnInternals=TRUE, plot=FALSE)
sd_train <- spec.mtm(d_train, Ftest=TRUE, returnInternals=TRUE, plot=FALSE)
sz_train <- spec.mtm(z_train, Ftest=TRUE, returnInternals=TRUE, plot=FALSE)

freq <- seq(1,nrow(sh_train$mtm$eigenCoefs),1) # frequencies at which to calculate the transfer function

tf <- get_tf_all(sh_train$mtm$eigenCoefs, sd_train$mtm$eigenCoefs, sz_train$mtm$eigenCoefs, s$mtm$eigenCoefs, freq)

#predict current ------------------
#spectral estimates of test data
sh_test <- spec.mtm(h_test, Ftest=TRUE, returnInternals=TRUE, plot=FALSE)
sd_test <- spec.mtm(d_test, Ftest=TRUE, returnInternals=TRUE, plot=FALSE)
sz_test <- spec.mtm(z_test, Ftest=TRUE, returnInternals=TRUE, plot=FALSE)


F_H <- matrix(ncol = length(freq), nrow = 7) # initialize empty matrix for sampling at given frequencies
F_D <- matrix(ncol = length(freq), nrow = 7) # initialize empty matrix for sampling at given frequencies
F_Z <- matrix(ncol = length(freq), nrow = 7) # initialize empty matrix for sampling at given frequencies


for(j in 1:length(freq)){  
  F_H[,j] <- t(sh_test$mtm$eigenCoefs[freq[j],]) #take row corresponding to jth chosen frequency
  F_D[,j] <- t(sd_test$mtm$eigenCoefs[freq[j],]) #take row corresponding to jth chosen frequency
  F_Z[,j] <- t(sz_test$mtm$eigenCoefs[freq[j],]) #take row corresponding to jth chosen frequency
}

pred_F_A <- matrix(ncol = length(freq), nrow = 7) # initialize empty matrix for sampling at given frequencies

for(i in 1:length(freq)){
  F_block <- as.matrix(data.frame(F_H[,i], F_D[,i], F_Z[,i]))  
  pred_F_A[,i] <- F_block %*% tf[i,]
}

prediction <- matrix(ncol = 7, nrow = nrow(sh_test$mtm$eigenCoefs))

for(j in 1:length(freq)){
  prediction[freq[j],] <- t(pred_F_A[,j])
}

#cmv's of prediction -----------
#complex mean values:
slep <- dpss(n = length(a_test), k = k, nw = 4, returnEigenvalues = FALSE)$v #slepian tapers

U_kzero <- mvfft(slep)[1, ] # You only want the zeroeth frequency

if (k >= 2){
  U_kzero[seq(2,k,2)] <- 0  # only want the even tapers... see P&W
}
ssqU_kzero <- sum(U_kzero^2)

all_pred_cmv <- (prediction %*% U_kzero) / ssqU_kzero

#F test on prediction ----------
k=7
pred_ftest <- rep(0, nrow(prediction)) #initialize vector of f test results

for(i in 1:nrow(prediction)){
  # formula from Thomson, 1982 (equation 13.10)
  denom = 0;

  for(j in 1:k){
    denom <- denom + abs(prediction[i,j] - all_pred_cmv[i]*U_kzero[j])^2
  }
  
  pred_ftest[i] <- Re(((k-1)*ssqU_kzero*abs(all_pred_cmv[i])^2)/denom)
    
}

pred_freqIndex <- prediction_fmax(pred_ftest, k, 0.9)

#reconstruct prediction----------

pred_cmv <- rep(0, nrow(prediction))
pred_cmv[pred_freqIndex] <- all_pred_cmv[pred_freqIndex] #only take frequencies from f test
pred_cmv[length(all_pred_cmv)-pred_freqIndex+2] <- Conj(all_pred_cmv[pred_freqIndex])

#for(i in 1:length(pred_cmv)){
#  pred_cmv[length(pred_cmv)-i] <- Conj(all_pred_cmv[i])  
#}

pred_recon <- Re(fft(pred_cmv, inverse=TRUE)[1:length(a_test)])

#reconstruction of actual data ---------------
# the fourier transform is complex, however since we started with real data, 
# the imaginary part will be zero, so just take the real part.
recon <- Re(fft(cmv, inverse=TRUE)[1:length(a2)])
fullRecon <- recon + cubeTrend

#plot with original data
par(mfrow=c(2,1))
plot(a2, type='l', col='grey', lwd=2, main="Reconstruction of Zero Mean Data")
lines(recon, col='blue')
plot(a_train, type='l', col='grey', lwd=2, main="Reconstruction of Original Data")
lines(fullRecon, col='blue')

#plot prediction------------

par(mfrow=c(2,1))
plot(a_train, type='l', col='grey', lwd=2, main="Reconstruction of Training Data")
lines(fullRecon, col='blue')
plot(a_test, type='l', col='grey', lwd=2, main="Reconstruction of Prediction on Testing Data")
lines(pred_recon+cubeTrend, col='red')

