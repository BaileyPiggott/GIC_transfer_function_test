# test method of creating a transfer function using simulated data

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

sim_data <- as.data.frame(read.table("sim_data.txt"))
colnames(sim_data) <- c('H', 'D', 'Z', 'A')

# create spectral estimates --------------------

N = 300 # number of data points per block

spec_H_block <- get_spec_mtm(sim_data$H, N) 
spec_D_block <- get_spec_mtm(sim_data$D, N)
spec_Z_block <- get_spec_mtm(sim_data$Z, N)
spec_A_block <- get_spec_mtm(sim_data$A, N)
# get_spec_mtm joins a frequency vector to the front of the data frame


# make transfer function ---------------------

freq <- seq(1,nrow(spec_H_block),2) # frequencies at which to calculate the transfer function

tf <- get_tf_ind(spec_H_block, spec_D_block, spec_Z_block, spec_A_block, freq)

tf_all <- get_tf_all(spec_H_block, spec_D_block, spec_Z_block, spec_A_block, freq)

#plot transfer functions ----------------------------

#individual plots assuming A has dependence on only given variable
ggplot(data = tf, aes(x = freq, y = val)) +
  facet_grid(tf~.) +
  geom_line() +
  labs(title = paste0("Transfer Functions\n Block Length = ", N, " samples\n", "Sampled at ", length(freq), " Frequencies"), x = "Frequency", y = "")

#facet plot of all three transfer functions
# assuming A is dependant on all variables

tf_three <- data.frame('freq' = freq, 'H' = tf_all[1,], 'D' = tf_all[2,], 'Z' = tf_all[3,]) %>% 
  gather(tf, val, 2:4)


ggplot(data = tf_three, aes(x = freq, y = val)) +
  facet_grid(tf~.) +
  geom_line() +
  labs(title = paste0("Transfer Functions\n Block Length = ", N, " samples\n", "Sampled at ", length(freq), " Frequencies"), x = "Frequency", y = "")
