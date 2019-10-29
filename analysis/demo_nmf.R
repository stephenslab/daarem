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
source("../code/nmf.R")

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# LOAD DATA
# ---------
cat("Loading droplet read counts.\n")
load("../data/droplet.RData")
droplet <- as.matrix(droplet)

# Randomly initialize the factorization.
n <- nrow(droplet)
m <- ncol(droplet)
A <- matrix(runif(n*k),n,k)
B <- matrix(runif(m*k),k,m)

# RUN BASIC EM UPDATES
# --------------------
cat("Running the multiplcative (EM) updates.\n")
fit1 <- betanmf(droplet,A,B,20)

# RUN ACCELERATED EM UPDATES
# --------------------------
cat("Running the accelerated EM updates.\n")
# TO DO.
