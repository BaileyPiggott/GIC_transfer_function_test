# method of creating a transfer function using simulated data
# 3 dimensional magnetic field example

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

#load data and functions and libraries------------
library(tidyr)
library(magrittr)
library(dplyr)
library(multitaper)
library(MASS) # need for pseudoinverse
library(ggplot2)

source('get_spec.R') #function to generate spectral estimate
source('get_spec_mtm.R') # generate spectral estimate using multitaper method
source('get_tf.R') # function to generate transfer function
source('get_tf_all.R')

#get data and separate into training and testing -----------------
sim_data <- as.data.frame(read.table("MA6.txt"))
sim_data <- sim_data[,-1] #MA4 and MA6 have an extra column
colnames(sim_data) <- c('H', 'D', 'Z', 'A')

##70% of the data for training; 30% for testing
#training_size <- ceiling(0.7 * nrow(sim_data))
#test_size <- nrow(sim_data)-training_size

#train_data <- sim_data[1:training_size,]
#test_data <- sim_data[(training_size+1):nrow(sim_data),]
train_data <- sim_data[1:1316,]
test_data <- sim_data[1317:2632,]

## plot simulated data --------------------------
# plot_sim <- sim_data %>% 
#   mutate(time = (1:nrow(sim_data))) %>% # add time column for plotting
#   gather(type, val, H:A) # convert to long form to plot in facets
# 
# ggplot(data = plot_sim, aes(x = time, y = val)) +
#   facet_grid(type~., scales = "free_y") +
#   geom_line() +
#   coord_cartesian(xlim = c(0,3000)) +
#   labs(title = "Simulated Data", x = "Time", y = "")+
#   theme(
#     title = element_text(size = 13),
#     axis.text = element_text(size  =10),
#     strip.text = element_text(size = 13, face = 'bold')
#   )

# create spectral estimates --------------------
block_N = 188 # number of data points per block
nw <- 4
k <- 7

# not using multitaper method
#spec_H_block <- get_spec(train_data$H, block_N)
#spec_D_block <- get_spec(train_data$D, block_N)
#spec_Z_block <- get_spec(train_data$Z, block_N)
#spec_A_block <- get_spec(train_data$A, block_N)

#using multitaper method
spec_H_block <- get_spec_mtm(train_data$H, nw, k, block_N) 
spec_D_block <- get_spec_mtm(train_data$D, nw, k, block_N)
spec_Z_block <- get_spec_mtm(train_data$Z, nw, k, block_N)
spec_A_block <- get_spec_mtm(train_data$A, nw, k, block_N)

# make transfer function ---------------------
freq <- seq(1,nrow(spec_H_block),1) # frequencies at which to calculate the transfer function

tf <- get_tf_all(spec_H_block, spec_D_block, spec_Z_block, spec_A_block, freq)

# plot transfer function --------------
# convert transfer function to matrix for graphing
p=1/2*nrow(tf)
mag <- abs(tf)
phase <- atan2(Im(tf), Re(tf))

tf_x <- data.frame(freq = seq(0,0.5,0.5/p), Magnitude = mag[1:(p+1),1], Phase = phase[1:(p+1),1]) %>%
  mutate(component = 'H')
tf_y <- data.frame(freq = seq(0,0.5,0.5/p), Magnitude = mag[1:(p+1),2], Phase = phase[1:(p+1),2]) %>%
  mutate(component = 'D')
tf_z <- data.frame(freq = seq(0,0.5,0.5/p), Magnitude = mag[1:(p+1),3], Phase = phase[1:(p+1),3]) %>%
  mutate(component = 'Z')

tf_plot <- bind_rows(tf_x, tf_y, tf_z) %>% 
  gather(type, val, Magnitude:Phase)

ggplot(data = tf_plot, aes(x = freq, y = val)) +
  facet_grid(type~component, scale = "free_y") +
  geom_line() +
  stat_smooth(method = "loess", formula = y ~ x, size = 0., se = "FALSE", colour = "red") +
  labs(title = paste0("Transfer Functions\nNumber of Blocks = ", (nrow(sim_data)-1)/(2*block_N)), x = "Frequency", y = "")+
  theme(
    title = element_text(size = 13),
    axis.text = element_text(size  =10),
    strip.text = element_text(size = 13, face = 'bold')
  )

#prediction ----------------------------

#using multitaper method:
test_H_block <- get_spec_mtm(test_data$H, nw, k, block_N) 
test_D_block <- get_spec_mtm(test_data$D, nw, k, block_N)
test_Z_block <- get_spec_mtm(test_data$Z, nw, k, block_N)
test_A_block <- get_spec_mtm(test_data$A, nw, k, block_N)


F_test_H <- matrix(ncol = length(freq), nrow = (ncol(test_H_block))) # initialize empty matrix for sampling at given frequencies
F_test_D <- matrix(ncol = length(freq), nrow = (ncol(test_D_block))) # initialize empty matrix for sampling at given frequencies
F_test_Z <- matrix(ncol = length(freq), nrow = (ncol(test_Z_block))) # initialize empty matrix for sampling at given frequencies


for(j in 1:length(freq)){  
  F_test_H[,j] <- t(test_H_block[freq[j],]) #take row corresponding to jth chosen frequency
  F_test_D[,j] <- t(test_D_block[freq[j],]) #take row corresponding to jth chosen frequency
  F_test_Z[,j] <- t(test_Z_block[freq[j],]) #take row corresponding to jth chosen frequency
}

F_predict_A <- matrix(ncol = length(freq), nrow = (ncol(test_A_block))) # initialize empty matrix for sampling at given frequencies

for(i in 1:length(freq)){
  F_block <- as.matrix(data.frame(F_test_H[,i], F_test_D[,i], F_test_Z[,i]))  
  F_predict_A[,i] <- F_block %*% tf[i,]
}

spec_pred_A <- matrix(ncol = (ncol(test_A_block)), nrow = (nrow(test_A_block)))

for(j in 1:length(freq)){
  spec_pred_A[freq[j],] <- t(F_predict_A[,j])
}

#convert prediction to time domain -----------

#complex mean values:
slep <- dpss(n = block_N, k = k, nw = nw, returnEigenvalues = FALSE)$v #slepian tapers

U_kzero <- mvfft(slep)[1, ] # You only want the zeroeth frequency

if (k >= 2){
  U_kzero[seq(2,k,2)] <- 0  # only want the even tapers... see P&W
}
ssqU_kzero <- sum(U_kzero^2)

#separate blocks and reconstruct them individually **********************
test <- spec_pred_A
spec_pred_A <- spec_pred_A[,1:k]

all_pred_cmv <- (spec_pred_A %*% U_kzero) / ssqU_kzero

#F test on prediction 

pred_ftest <- rep(0, nrow(spec_pred_A)) #initialize vector of f test results

for(i in 1:nrow(spec_pred_A)){
  # formula from Thomson, 1982 (equation 13.10)
  denom = 0;
  
  for(j in 1:k){
    denom <- denom + abs(spec_pred_A[i,j] - all_pred_cmv[i]*U_kzero[j])^2
  }
  pred_ftest[i] <- Re(((k-1)*ssqU_kzero*abs(all_pred_cmv[i])^2)/denom)
}

pred_freqIndex <- prediction_fmax(pred_ftest, k, 0.9)

#reconstruct prediction----------

pred_cmv <- rep(0, nrow(spec_pred_A))
pred_cmv[pred_freqIndex] <- all_pred_cmv[pred_freqIndex] #only take frequencies from f test
pred_cmv[length(all_pred_cmv)-pred_freqIndex+2] <- Conj(all_pred_cmv[pred_freqIndex])

pred_recon <- Re(fft(pred_cmv, inverse=TRUE)[1:block_N])

#plot ---------------

plot(test_data$A[1:block_N], type='l', col = "grey", lwd=2,  main= "Prediction vs. Measured Data")
lines(pred_recon, col='red')


