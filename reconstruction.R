## recontruction of data


# removing a cubic trend (want this for the reconstruction part)
sim_data <- as.data.frame(read.table("MA6.txt"))
d <- sim_data$V5 #current column

x <- 1:length(d)
fit <- lm(d ~ x + I(x^2) + I(x^3))
cubeTrend <- fitted.values(fit)
d2 <- d - cubeTrend



# plot-check----------------
## top: plot with cubic fit overlay
## bottom: residuals after removing cubic fit.
par(mfrow=c(2,1), mar=c(4,4,1,1), mgp=c(2.5,1,0))
plot(d, type='l', main = "Original Data with Cubic Trend")
lines(cubeTrend, col='red', lwd=2)
plot(d2, type='l', ylab="Residuals", main = "Residuals") # want 0-mean for reconstruction


#calculate spectrum ----------------------
require('multitaper')

# calculate the spectrum (returnInternals gives us the complex-mean-values (cmv's))
# cmv's are an average of the eigencoefficients (windowed-fourier-transformed data)
# if you average the squared eigencoefficients, you get the spectrum.

s <- spec.mtm(d2, Ftest=TRUE, returnInternals=TRUE, plot=FALSE)

freqIdx <- findLocalFMax(s, 0.9)# find frequencies of large F-values ("most sinusoidal" frequencies)

# set up the CMV array (spec.mtm() only returns the positive frequencies due to symmetry)
cmv <- rep(0, s$mtm$nFFT)
cmv[freqIdx] <- s$mtm$cmv[freqIdx]
cmv[s$mtm$nFFT-freqIdx+2] <- Conj(s$mtm$cmv[freqIdx])

# the "negative" side of the cmv's actually go in the right side of the array and must be
# conjugated (due to the minus sign on the frequency)

#reconstruction ---------------
recon <- Re(fft(cmv, inverse=TRUE)[1:length(d2)])
# the fourier transform is complex, however since we started with real data, 
# the imaginary part will be zero, so just take the real part.

fullRecon <- recon + cubeTrend

par(mfrow=c(2,1))

plot(d2, type='l', col='grey', lwd=2, main="Reconstruction of Zero Mean Data")
lines(recon, col='blue')
plot(d, type='l', col='grey', lwd=2, main="Reconstruction of Original Data")
lines(fullRecon, col='blue')


#Function -------------------
# Finds local maxes in the Ftest given a spec.mtm object.
# based on a probability cutoff

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
