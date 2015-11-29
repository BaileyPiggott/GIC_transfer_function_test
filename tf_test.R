# method of creating a transfer function using simulated data

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

# plot simulated data --------------------------
plot_sim <- sim_data %>% 
  mutate(time = (1:nrow(sim_data))) %>% # add time column for plotting
  gather(type, val, H:A) # convert to long form to plot in facets

ggplot(data = plot_sim, aes(x = time, y = val)) +
  facet_grid(type~.) +
  geom_line() +
  labs(title = "Simulated Data", x = "Time", y = "")

# create spectral estimates --------------------

block_N = 200 # number of data points per block
nw <- 4
k <- 7

spec_H_block <- get_spec_mtm(sim_data$H, nw, k, block_N) 
spec_D_block <- get_spec_mtm(sim_data$D, nw, k, block_N)
spec_Z_block <- get_spec_mtm(sim_data$Z, nw, k, block_N)
spec_A_block <- get_spec_mtm(sim_data$A, nw, k, block_N)

# make transfer function ---------------------

freq <- seq(1,nrow(spec_H_block),2) # frequencies at which to calculate the transfer function

tf_H <- get_tf_ind(spec_H_block, spec_A_block, freq)
tf_D <- get_tf_ind(spec_D_block, spec_A_block, freq)
tf_Z <- get_tf_ind(spec_Z_block, spec_A_block, freq)

#tf_all <- get_tf_all(spec_H_block, spec_D_block, spec_Z_block, spec_A_block, freq)

# plot transfer functions --------------

# independent transfer functions
ggplot(data = tf_H, aes(x = freq, y = val)) +
  facet_grid(type~., scale = "free_y") +
  geom_line() +
  stat_smooth(method = "loess", formula = y ~ x, size = 1, se = "FALSE", colour = "red") +
  labs(title = paste0("Transfer Function\nNumber of Blocks = ", nrow(sim_data)/block_N), x = "Frequency", y = "")


#-----------------------

#facet plot of all three transfer functions
# assuming A is dependant on all variables

tf_three <- data.frame('freq' = freq, 'H' = tf_all[1,], 'D' = tf_all[2,], 'Z' = tf_all[3,]) %>% 
  gather(tf, val, 2:4)


ggplot(data = tf_three, aes(x = freq, y = val)) +
  facet_grid(tf~.) +
  geom_line() +
  labs(title = paste0("Transfer Functions\n Block Length = ", N, " samples\n", "Sampled at ", length(freq), " Frequencies"), x = "Frequency", y = "")
