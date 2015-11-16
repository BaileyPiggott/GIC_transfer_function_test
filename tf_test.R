# test method of creating a transfer function using simulated data
#load data and libraries------------
library(tidyr)
library(magrittr)
library(dplyr)
library(multitaper)

data <- read.table("sim_data.txt")
colnames(data) <- c('H', 'D', 'Z', 'A')

# create periodogram --------------------
#plot(data$A, type = 'l')

N = 200
H <- data$H[1:N]
A <- data$A[1:N]

spec_H <- (1/N)*abs(fft(H))^2 # spectral estimate
spec_A <- (1/N)*abs(fft(A))^2 # spectral estimate

# make transfer function ---------------------




