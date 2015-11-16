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
N = 200 # number of data points per block
#initialize empty matrices where each block is a column
H_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))
D_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))
Z_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))
A_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))

spec_H_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))
spec_D_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))
spec_Z_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))
spec_A_block <- data.frame(matrix(ncol = (nrow(data)/N), nrow = N))

# create blocks of the data
for(j in 1:1){#(nrow(data)/N)){
  H_block[,j]<- data$H[(1+(j-1)*N):(j*N)]
  D_block[,j]<- data$D[(1+(j-1)*N):(j*N)]
  Z_block[,j]<- data$Z[(1+(j-1)*N):(j*N)]
  A_block[,j]<- data$A[(1+(j-1)*N):(j*N)]  

#spectral estimate of each block
 spec_H_block[,j] <- (1/N)*abs(fft(H_block[,j]))^2
 spec_D_block[,j] <- (1/N)*abs(fft(D_block[,j]))^2
 spec_Z_block[,j] <- (1/N)*abs(fft(Z_block[,j]))^2
 spec_A_block[,j] <- (1/N)*abs(fft(A_block[,j]))^2
}


# make transfer function ---------------------

freq <- c(1:N)
M=vector(mode = 'numeric', length = length(freq))

for(i in 1:length(freq)){
F_A <- as.matrix(data.frame('F' = c(spec_A_1[freq[i]], spec_A_2[freq[i]], spec_A_3[freq[i]])))
F_H <- as.matrix(data.frame('F' = c(spec_H_1[freq[i]], spec_H_2[freq[i]], spec_H_3[freq[i]])))

M[i] <- ginv(F_H) %*% F_A
}

plot(M, type = 'l')

