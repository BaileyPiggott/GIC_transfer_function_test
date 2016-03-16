# calculate and plot coherence between the current and magnetic field components

#load functions ----------
library(tidyr)
library(magrittr)
library(dplyr)
library(multitaper)
library(MASS) # need for pseudoinverse
library(ggplot2)
library(grid)

select <- dplyr::select # conflicts with 'select' in MASS

# load data ------
geomag <- readRDS("boulder_1999-2012_15min_final.rds")
inducedCurrent <- readRDS("inducedCurrent-1999-2012-15min.rds")

data <- inner_join(geomag, inducedCurrent, by = "date") %>% select(X,Y,Z,value)
colnames(data) <- c("H", "D", "Z", "A")

nw = 4
k = 7

#calculate coherence----------

H <- spec.mtm(data[1:1000,]$H, nw=nw, k=k, returnInternals = TRUE, plot= FALSE)
D <- spec.mtm(data[1:1000,]$D, nw=nw, k=k, returnInternals = TRUE, plot= FALSE)
Z <- spec.mtm(data[1:1000,]$Z, nw=nw, k=k, returnInternals = TRUE, plot= FALSE)
A <- spec.mtm(data[1:1000,]$A, nw=nw, k=k, returnInternals = TRUE, plot= FALSE)

Hcoh <- mtm.coh(H, A, plot = TRUE)
Dcoh <- mtm.coh(D, A, plot = TRUE)
Zcoh <- mtm.coh(Z, A, plot = TRUE)

# convert into one data frame for facet plots
coherence <- data.frame(Hcoh$freq,  Hcoh$msc, Dcoh$msc, Zcoh$msc)
colnames(coherence) <- c("freq", "H", "D", "Z")

coherence <- coherence %>% gather("type", "MSC", H:Z)

# plot -------------------
ggplot(data = coherence,aes(x= freq, y = MSC)) +
  facet_grid(.~type) +
  geom_line()+
  coord_cartesian( ylim = c(0, 1)) + 
  theme(
    axis.line = element_line("grey"), 
    panel.grid.major = element_line("grey"),
    panel.background = element_rect("white"),
    axis.ticks = element_blank(), # remove x axis ticks
    strip.text = element_text(size = 13, face = 'bold'),
    panel.margin = unit(1.25, "lines"), # space between facets
    panel.border = element_rect(fill=NA, "grey"),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12) #size of x axis labels
  ) +
  labs(title = "Coherence between Magnetic Field and Measured Current", x = "Frequency", y = "Magnitude Squared Coherence")
