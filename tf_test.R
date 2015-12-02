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
  coord_cartesian(xlim = c(0,3000)) +
  labs(title = "Simulated Data", x = "Time", y = "")+
  theme(
    title = element_text(size = 13),
    axis.text = element_text(size  =10),
    strip.text = element_text(size = 13, face = 'bold')
  )

# create spectral estimates --------------------

block_N = 3000 # number of data points per block
nw <- 4
k <- 7

# not using multitaper method
spec_H_block <- get_spec(sim_data$H, block_N)
spec_D_block <- get_spec(sim_data$D, block_N)
spec_Z_block <- get_spec(sim_data$Z, block_N)
spec_A_block <- get_spec(sim_data$A, block_N)

#using multitaper method
spec_H_block <- get_spec_mtm(sim_data$H, nw, k, block_N) 
spec_D_block <- get_spec_mtm(sim_data$D, nw, k, block_N)
spec_Z_block <- get_spec_mtm(sim_data$Z, nw, k, block_N)
spec_A_block <- get_spec_mtm(sim_data$A, nw, k, block_N)

# make transfer function ---------------------

freq <- seq(1,nrow(spec_H_block),2) # frequencies at which to calculate the transfer function

tf_all <- get_tf_all(spec_H_block, spec_D_block, spec_Z_block, spec_A_block, freq)

# plot transfer functions --------------

ggplot(data = tf_all, aes(x = freq, y = val)) +
  facet_grid(type~component, scale = "free_y") +
  geom_line() +
  stat_smooth(method = "loess", formula = y ~ x, size = 0., se = "FALSE", colour = "red") +
  labs(title = paste0("Transfer Functions\nNumber of Blocks = ", nrow(sim_data)/block_N), x = "Frequency", y = "")+
  theme(
    title = element_text(size = 13),
    axis.text = element_text(size  =10),
    strip.text = element_text(size = 13, face = 'bold')
  )

