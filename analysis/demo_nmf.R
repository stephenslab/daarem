# TO DO: Explain what this script is for, and how to use it.
#
# X = A*B
#

# SCRIPT PARAMETERS
# -----------------
# The rank of the matrix factorization (number of "topics).
k <- 3

# SET UP ENVIRONMENT
# ------------------
# Load some packages and function definitions used in the code below.
library(Matrix)
library(daarem)
source("../code/nmf.R")

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# LOAD DATA
# ---------
# cat("Loading droplet read counts.\n")
# load("../data/droplet.RData")
# X <- as.matrix(droplet)
# rm(droplet)
n <- 300
m <- 400
F <- matrix(runif(m*k)/3,m,k)
L <- matrix(runif(n*k)/3,n,k)
X <- matrix(rpois(n*m,tcrossprod(L,F)),n,m)

# Randomly initialize the factorization.
# n <- nrow(X)
# m <- ncol(X)
A <- matrix(runif(n*k),n,k)
B <- matrix(runif(m*k),k,m)

# RUN BASIC EM UPDATES
# --------------------
cat("Running the multiplcative (EM) updates.\n")
fit1 <- betanmf(X,A,B,100)

# RUN ACCELERATED EM UPDATES
# --------------------------
cat("Running the accelerated EM updates.\n")
fit2 <- daarbetanmf(X,A,B,100)

library(ggplot2)
library(cowplot)
pdat <-
  rbind(data.frame(iter = 1:100,objective = fit1$value,method = "EM"),
        data.frame(iter = 1:100,objective = fit2$value,method = "daarem"))
ggplot(pdat,aes(x = iter,y = objective,color = method,linetype = method)) +
  geom_line(size = 1)
