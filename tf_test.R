# method of creating a transfer function using simulated data
# 3 dimensional magnetic field example

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

##get data and separate into training and testing -----------------
#data <- as.data.frame(read.table("MA6.txt"))
#data <- data[,-1] #MA4 and MA6 have an extra column
#colnames(data) <- c('H', 'D', 'Z', 'A')


# geomag <- readRDS("boulder_1999-2012_15min_final.rds")
# inducedCurrent <- readRDS("inducedCurrent-1999-2012-15min.rds")
# 
# data <- inner_join(geomag, inducedCurrent, by = "date") %>% select(X,Y,Z,value)
# colnames(data) <- c("H", "D", "Z", "A")
# 
# 
# ##70% of the data for training; 30% for testing
# training_size <- ceiling(0.7 * nrow(data))
# 
# train_data <- data[1:training_size,]
# test_data <- data[(training_size+1):(nrow(data)),]


#fit <- lm(train_data ~ x + I(x^2) + I(x^3))
#cubeTrend <- fitted.values(fit)
#a2 <- train_data - cubeTrend

# test bad day data -------------

geomag <- readRDS("boulder_1999-2012_15min_final.rds")
inducedCurrent <- readRDS("inducedCurrent-1999-2012-15min.rds")

data <- inner_join(geomag, inducedCurrent, by = "date") %>% select(X,Y,Z,value)
colnames(data) <- c("H", "D", "Z", "A")


maxInd <- which(abs(data$A) > 1)
bad_day <- as.data.frame(matrix(nrow = length(maxInd), ncol = 4))


for(i in 1:length(maxInd)){
  
  bad_day[i, ] <- data[maxInd[i], ]
  
}

colnames(bad_day) <- colnames(data)
train_data <- bad_day

test_data <- data[2001:3000, ]

## plot data --------------------------
# plot_sim <- data %>% 
#   mutate(time = (1:nrow(data))) %>% # add time column for plotting
#   gather(type, val, H:A) # convert to long form to plot in facets
# 
# ggplot(data = plot_sim, aes(x = time, y = val)) +
#   facet_grid(type~., scales = "free_y") +
#   geom_line() +
#   coord_cartesian(xlim = c(0,3000)) +
#   labs(title = "Simulated Data", x = "Time", y = "")+
#   theme(
#     title = element_text(size = 13),
#     axis.text = element_text(size  =10),
#     strip.text = element_text(size = 13, face = 'bold')
#   )

# create spectral estimates --------------------
block_N = 500 # number of data points per block; 
# *****N MUST be a factor of number of samples in test data***
nw <- 4
k <- 7

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
time_data$index <- as.numeric(rownames(time_data))

time_data <- time_data %>% gather("type", "a", measured:predicted)

#plot ---------------

ggplot(data = time_data,aes(x= index, y = a, colour = type, size = type)) +
  geom_line()+
  #coord_cartesian( ylim = c(-2.5, 2.5)) + 
  theme(
    axis.line = element_line("grey"), 
    panel.grid.major.y = element_line("grey"),
    panel.grid.major.x = element_blank(), # remove vertical white lines
    panel.background = element_rect("white"),
    axis.ticks.x = element_blank(), # remove x axis ticks
    plot.title = element_text(size = 15),
    legend.key = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12) #size of x axis labels
  ) +
  labs(title = "Measured Current vs. Predicted", x = "Time (15 min interval)", y = "Induced Current (A)") +
  scale_colour_manual(
    values = c("grey", "red"), 
    name = "" ) +
  scale_size_manual(
    values = c(1.2,0.7), 
    name = "")
  
  
