# test method of creating a transfer function using simulated data

#load data and functions and libraries------------
library(tidyr)
library(magrittr)
library(dplyr)
library(multitaper)
library(MASS) # need for pseudoinverse
library(ggplot2)

source('get_spec.R') #function to generate spectral estimate
source('get_tf.R') # function to generate transfer function

sim_data <- as.data.frame(read.table("sim_data.txt"))
colnames(sim_data) <- c('H', 'D', 'Z', 'A')

# create spectral estimates --------------------

N = 300 # number of data points per block

spec_H_block <- get_spec(sim_data$H, N)
spec_D_block <- get_spec(sim_data$D, N)
spec_Z_block <- get_spec(sim_data$Z, N)
spec_A_block <- get_spec(sim_data$A, N)

# make transfer function ---------------------

freq <- seq(1,(N/2+1),2) # frequencies at which to calculate the transfer function

tf <- get_tf_ind(spec_H_block, spec_D_block, spec_Z_block, spec_A_block, freq)

#tf_all <- get_tf_all()

#plot transfer functions ----------------------------

#individual plots assuming A has dependence on only given variable
ggplot(data = tf, aes(x = freq, y = val)) +
  facet_grid(tf~.) +
  geom_line() +
  labs(title = paste0("Transfer Functions\n Block Length = ", N, " samples\n", "Sampled at ", length(freq), " Frequencies"), x = "Frequency", y = "")


#facet plot of all three transfer functions
# assuming A is dependant on all variables

tf_all <- data.frame('freq' = freq, 'H' = M_test[,1], 'D' = M_test[,2], 'Z' = M_test[,3]) %>% 
  gather(tf, val, 2:4)

ggplot(data = tf_all, aes(x = freq, y = val)) +
  facet_grid(tf~.) +
  geom_line() +
  labs(title = paste0("Transfer Functions\n Block Length = ", N, " samples"), x = "Frequency", y = "")

