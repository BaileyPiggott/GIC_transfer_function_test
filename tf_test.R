# test method of creating a transfer function using simulated data
#load data and libraries------------
library(tidyr)
library(magrittr)
library(dplyr)
library(multitaper)
library(MASS) # need for pseudoinverse

source("get_tf.R") # function to get transfer function at given freq

data <- read.table("sim_data.txt")
colnames(data) <- c('H', 'D', 'Z', 'A')

# create periodogram --------------------
#plot(data$A, type = 'l')

N = 200
H_1 <- data$H[1:N]
A_1 <- data$A[1:N]

H_2 <- data$H[(N+1):(2*N)]
A_2 <- data$A[(N+1):(2*N)]

H_3 <- data$H[(2*N+1):(3*N)]
A_3 <- data$A[(2*N+1):(3*N)]

spec_H_1 <- (1/N)*abs(fft(H_1))^2 # spectral estimate
spec_A_1 <- (1/N)*abs(fft(A_1))^2 # spectral estimate

spec_H_2 <- (1/N)*abs(fft(H_2))^2 # spectral estimate
spec_A_2 <- (1/N)*abs(fft(A_2))^2 # spectral estimate

spec_H_3 <- (1/N)*abs(fft(H_3))^2 # spectral estimate
spec_A_3 <- (1/N)*abs(fft(A_3))^2 # spectral estimate

# make transfer function ---------------------

freq <- c(1:N)
M=vector(mode = 'numeric', length = length(freq))

for(i in 1:length(freq)){
F_A <- as.matrix(data.frame('F' = c(spec_A_1[freq[i]], spec_A_2[freq[i]], spec_A_3[freq[i]])))
F_H <- as.matrix(data.frame('F' = c(spec_H_1[freq[i]], spec_H_2[freq[i]], spec_H_3[freq[i]])))

M[i] <- ginv(F_H) %*% F_A
}

plot(M, type = 'l')

