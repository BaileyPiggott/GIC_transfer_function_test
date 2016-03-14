
#load functions and libraries------------
library(tidyr)
library(magrittr)
library(dplyr)
library(multitaper)
library(MASS) # need for pseudoinverse
library(ggplot2)

select <- dplyr::select # conflicts with 'select' in MASS

source('get_spec.R') #function to generate spectral estimate
source('get_spec_mtm.R') # generate spectral estimate using multitaper method
source('get_tf.R') # function to generate transfer function
source('get_tf_all.R')
source('recon_function.R')

# parameters to adjust ----------
block_N = 500 # number of data points per block; 
# *****N MUST be a factor of number of samples in test data***

current_threshold <- 1 #for "bad day" data
test_data_low <- 3001
test_data_high <- 3500 #interval of testing data

nw <- 4
k <- 7

#get data and separate into training and testing -----------------

geomag <- readRDS("boulder_1999-2012_15min_final.rds")
inducedCurrent <- readRDS("inducedCurrent-1999-2012-15min.rds")

data <- inner_join(geomag, inducedCurrent, by = "date") %>% select(X,Y,Z,value)
colnames(data) <- c("H", "D", "Z", "A")

bad_day <- as.data.frame(matrix(nrow = length(data), ncol = 4)) #initalize "bad day data"

j=1;
for (i in 1:nrow(data)){
  
  if(data$A[i] > current_threshold){
    bad_day[j, ] = data[i-1, ]
    j=j+1
    bad_day[j,] = data[i, ]
    j=j+1
    bad_day[j,] = data[i+1, ]
    j=j+1
  }
}
colnames(bad_day) <- colnames(data)

#choose training data:
train_data <- bad_day #uncomment to use bad day data as training data
#train_data <-data[4000:14068, ] #uncomment to use random data block as training data

test_data <- data[test_data_low:test_data_high, ] #testing data frame based on chosen thresholds

# # plot time domain data --------------------------
# plot_sim <- data[11150:11520, ] %>% 
#   mutate(time = (11150:11520)) %>% # add time column for plotting
#   gather(type, val, H:A) # convert to long form to plot in facets
# 
# plot_sim$time <- plot_sim$time*15/60/24
# 
# ggplot(data = plot_sim, aes(x = time, y = val)) +
#   facet_grid(type~., scales = "free_y") +
#   geom_line() +
#   labs(title = "Subset of Provided Data", x = "Time (Days)", y = "")+
#   theme(
#     plot.title = element_text(size = 16, face = "bold"),
#     axis.title = element_text(size = 14),
#     axis.text = element_text(size = 12), #size of x axis labels
#     strip.text = element_text(size = 13, face = 'bold')
#   )

# create spectral estimates --------------------

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

## plot transfer function --------------
# # convert transfer function to matrix for graphing
# p=1/2*nrow(tf)
# mag <- abs(tf)
# phase <- atan2(Im(tf), Re(tf))
# 
# tf_x <- data.frame(freq = seq(0,0.5,0.5/p), Magnitude = mag[1:(p+1),1], Phase = phase[1:(p+1),1]) %>%
#   mutate(component = 'H')
# tf_y <- data.frame(freq = seq(0,0.5,0.5/p), Magnitude = mag[1:(p+1),2], Phase = phase[1:(p+1),2]) %>%
#   mutate(component = 'D')
# tf_z <- data.frame(freq = seq(0,0.5,0.5/p), Magnitude = mag[1:(p+1),3], Phase = phase[1:(p+1),3]) %>%
#   mutate(component = 'Z')
# 
# tf_plot <- bind_rows(tf_x, tf_y, tf_z) %>% 
#   gather(type, val, Magnitude:Phase)
# 
# ggplot(data = tf_plot, aes(x = freq, y = val)) +
#   facet_grid(type~component, scale = "free_y") +
#   geom_line() +
#   stat_smooth(method = "loess", formula = y ~ x, size = 0., se = "FALSE", colour = "red") +
#   labs(title = paste0("Transfer Functions\nNumber of Blocks = ", (nrow(data)-1)/(2*block_N)), x = "Frequency", y = "")+
#   theme(
#     title = element_text(size = 13),
#     axis.text = element_text(size  =10),
#     strip.text = element_text(size = 13, face = 'bold')
#   )

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

#reconstruct prediction and format for plotting --------------------

prediction<- data.frame(0) # initialize prediction

for(i in 1:(nrow(test_data)/block_N)){
  
  pred_recon <- reconstruct(spec_pred_A[,((i-1)*k+1):(i*k)], block_N, k, nw)
  
  prediction<- bind_rows(prediction, pred_recon)
}
prediction <- prediction[-1, 2] #take out initialize row 

measured <- data.frame(test_data$A)

time_data <- bind_cols(measured, prediction)
colnames(time_data) <- c('measured', 'predicted')
time_data$index <- as.numeric(rownames(time_data))*15/60/24 #data is in 15 minute intervals

time_data <- time_data %>% gather("type", "a", measured:predicted)
time_data <- as.data.frame(time_data)

#plot ---------------

plot_title = paste0("Measured Current vs. Predicted")

ggplot(data = time_data,aes(x= index, y = a, colour = type, size = type)) +
  geom_line()+
  coord_cartesian( ylim = c(-0.7, 0.7)) + 
  theme(
    axis.line = element_line("grey"), 
    panel.grid.major.y = element_line("grey"),
    panel.grid.major.x = element_blank(), # remove vertical white lines
    panel.background = element_rect("white"),
    axis.ticks.x = element_blank(), # remove x axis ticks
    panel.border = element_rect(fill=NA, "grey"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.key = element_blank(),
    legend.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12) #size of x axis labels
  ) +
  labs(title = plot_title, x = "Time (Days)", y = "Induced Current (A)") +
  scale_colour_manual(
    values = c("grey70", "red"), 
    name = "" ) +
  scale_size_manual(
    values = c(1,0.7), 
    name = "")
  
  

#coherence and coherence plot -----------

real_data <- spec.mtm(test_data$A, nw=nw, k=k, returnInternals = TRUE, plot= FALSE)
prediction_data <- spec.mtm(prediction, nw=nw, k=k, returnInternals = TRUE, plot=FALSE)

resCoh <-mtm.coh(real_data, prediction_data)
plot(resCoh)
plot(resCoh$freq, resCoh$msc, type = "l")

coherence <- data.frame(resCoh$freq,  resCoh$msc)
colnames(coherence) <- c("freq", "msc")


ggplot(data = coherence,aes(x= freq, y = msc)) +
  geom_line()+
  coord_cartesian( ylim = c(0, 1)) + 
  theme(
    axis.line = element_line("grey"), 
    panel.grid.major.y = element_line("grey"),
    panel.grid.major.x = element_blank(), # remove vertical white lines
    panel.background = element_rect("white"),
    axis.ticks = element_blank(), # remove x axis ticks
    #panel.border = element_rect(fill=NA, "grey"),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12) #size of x axis labels
  ) +
  labs(title = "Coherence between Measured and Predicted Current", x = "Frequency", y = "Mean Squared Coherence")
