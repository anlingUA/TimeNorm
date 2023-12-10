# TimeNorm
# Normalization for time-series microbial count data

# group 2 is control group

# input raw data: example_data.csv
# 2 groups, 6 time points, 10 samples
# input data in wide format:
#      Group 1                    Group 2
# T1          T2   ...   Tk     T1          T2  ...    Tk   
# S1 S2 ..Sn  S1 S2 ..Sn        S1 S2 ..Sn  S1 S2 ..Sn

library(protoclust)
library(TANOVA)
library(data.table)
library(gplots)
library(RColorBrewer)


source("raida.R")
source("TimeNorm.R")
set.seed(2023)

TN <- function(wide.data)
{
  test1tmm <- step1Norm(wide.data,t=t,rep=rep,g=g)
  test2tmm <- step2Norm(test1tmm)
}


input_data <- read.csv("example_data.csv")
### define values
g <- 2  # group
t <- 6  # timepoint
rep <- 10 # samples


normaldata <- TN(input_data)
  #output normalized data from TN function
g1 <- do.call(cbind, normaldata$s2.c1.normal.dat)
g2 <- do.call(cbind, normaldata$s2.c2.normal.dat)
output <- cbind(g1,g2)

rg <- colorRampPalette(brewer.pal(8, "Reds"))(18)

png("Before_normalization.png")
heatmap.2(as.matrix(sqrt(input_data)), scale = "none", trace = "none",
          Rowv = F, Colv = F, breaks =seq(0,18,1), 
          dendrogram="none", col=rg,labRow = NA, labCol = NA)
dev.off()

png("After_normalization.png")
heatmap.2(as.matrix(sqrt(output)), scale = "none", trace = "none",
          Rowv = F, Colv = F, breaks =seq(0,18,1), 
          dendrogram="none", col=rg,labRow = NA, labCol = NA)
dev.off()
