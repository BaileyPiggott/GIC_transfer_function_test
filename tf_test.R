# test method of creating a transfer function using simulated data

#load data and libraries------------
library(tidyr)
library(magrittr)
library(dplyr)
library(multitaper)
library(MASS) # need for pseudoinverse

data <- read.table("sim_data.txt")
colnames(data) <- c('H', 'D', 'Z', 'A')

# create periodogram --------------------

# create multiple blocks of the data

N = 300 # number of data points per block

#initialize empty matrices where each block is a column
H_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))
D_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))
Z_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))
A_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))
# initialize empty matrices where spectral estimation of each block is a column
spec_H_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))
spec_D_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))
spec_Z_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))
spec_A_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))
# no windows used in spectral analysis yet

for(j in 1:(nrow(data)/N)){
#separate into blocks:  
  H_block[,j]<- data$H[(1+(j-1)*N):(j*N)]
  D_block[,j]<- data$D[(1+(j-1)*N):(j*N)]
  Z_block[,j]<- data$Z[(1+(j-1)*N):(j*N)]
  A_block[,j]<- data$A[(1+(j-1)*N):(j*N)]  

#spectral estimate of each block:
 spec_H_block[,j] <- (1/N)*abs(fft(H_block[,j]))^2
 spec_D_block[,j] <- (1/N)*abs(fft(D_block[,j]))^2
 spec_Z_block[,j] <- (1/N)*abs(fft(Z_block[,j]))^2
 spec_A_block[,j] <- (1/N)*abs(fft(A_block[,j]))^2
}

spec_H_block <- spec_H_block[1:(N/2+1),] #spectra are symmetric; only need half
spec_D_block <- spec_D_block[1:(N/2+1),]
spec_Z_block <- spec_Z_block[1:(N/2+1),]
spec_A_block <- spec_A_block[1:(N/2+1),]

# make transfer function ---------------------

freq <- c(1:(N/2+1)) # frequencies at which to calculate the transfer function

M_H=vector(mode = 'numeric', length = length(freq)) # empty matrix to store tf
M_D=vector(mode = 'numeric', length = length(freq)) # empty matrix to store tf
M_Z=vector(mode = 'numeric', length = length(freq)) # empty matrix to store tf

# rows correspond to blocks, and columns correspond to one frequency
F_H <- t(as.matrix(spec_H_block)) # sampling at every frequency
F_D <- t(as.matrix(spec_D_block))
F_Z <- t(as.matrix(spec_Z_block))
F_A <- t(as.matrix(spec_A_block))

#calculate and plot transfer function for each component (M)
for(i in 1:length(freq)){
  
  M_H[i] <- ginv(F_H[,i]) %*% F_A[,i]
  M_D[i] <- ginv(F_D[,i]) %*% F_A[,i]
  M_Z[i] <- ginv(F_Z[,i]) %*% F_A[,i]
}

plot(M_H, type = 'l')
plot(M_D, type = 'l')
plot(M_Z, type = 'l')

