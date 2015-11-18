# test method of creating a transfer function using simulated data

#load data and functions and libraries------------
library(tidyr)
library(magrittr)
library(dplyr)
library(multitaper)
library(MASS) # need for pseudoinverse
library(ggplot2)

source('get_spec.R') #function to generate spectral estimate

sim_data <- as.data.frame(read.table("sim_data.txt"))
colnames(sim_data) <- c('H', 'D', 'Z', 'A')

# create spectral estimates --------------------

N = 300 # number of data points per block

spec_H_block <- get_tf(sim_data$H, N)
spec_D_block <- get_tf(sim_data$D, N)
spec_Z_block <- get_tf(sim_data$Z, N)
spec_A_block <- get_tf(sim_data$A, N)


# make transfer function ---------------------

freq <- c(1:(N/2+1)) # frequencies at which to calculate the transfer function

M_H=vector(mode = 'numeric', length = length(freq)) # empty matrix to store tf
M_D=vector(mode = 'numeric', length = length(freq)) # empty matrix to store tf
M_Z=vector(mode = 'numeric', length = length(freq)) # empty matrix to store tf

M_test= matrix(ncol = 3, nrow = length(freq)) # empty matrix to store tf

# rows correspond to blocks, and columns correspond to one frequency
F_H <- t(as.matrix(spec_H_block)) # sampling at every frequency
F_D <- t(as.matrix(spec_D_block))
F_Z <- t(as.matrix(spec_Z_block))
F_A <- t(as.matrix(spec_A_block))



#calculate and plot transfer function for each component (M)
for(i in 1:length(freq)){

  F_block <- as.matrix(data.frame(F_H[,i], F_D[,i], F_Z[,i]))
  
  M_test[i,] <- ginv(F_block) %*% F_A[,i]
  
  M_H[i] <- ginv(F_H[,i]) %*% F_A[,i]
  M_D[i] <- ginv(F_D[,i]) %*% F_A[,i]
  M_Z[i] <- ginv(F_Z[,i]) %*% F_A[,i]
}


#plot transfer functions --------------

#combine transfer function data into one data frame with frequency
tf <- data.frame('freq' = freq, 'H' = M_H, 'D' = M_D, 'Z' = M_Z) %>% 
  gather(tf, val, 2:4)

#facet plot of all three transfer functions
#individual plots assuming A has dependence on only given variable
ggplot(data = tf, aes(x = freq, y = val)) +
  facet_grid(tf~.) +
  geom_line() +
  labs(title = paste0("Transfer Functions\n Block Length = ", N, " samples"), x = "Frequency", y = "")

  
tf_all <- data.frame('freq' = freq, 'H' = M_test[,1], 'D' = M_test[,2], 'Z' = M_test[,3]) %>% 
  gather(tf, val, 2:4)

#facet plot of all three transfer functions
# assuming A is dependant on all variables
ggplot(data = tf_all, aes(x = freq, y = val)) +
  facet_grid(tf~.) +
  geom_line() +
  labs(title = paste0("Transfer Functions\n Block Length = ", N, " samples"), x = "Frequency", y = "")

