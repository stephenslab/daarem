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
library(NNLM)
library(fastTopics)
source("../code/nmf.R")

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# SIMULATE DATA
# -------------
n <- 300
m <- 400
a <- 2 # 1/3
F <- matrix(a*runif(m*k),m,k)
L <- matrix(a*runif(n*k),n,k)
X <- matrix(rpois(n*m,tcrossprod(L,F)),n,m)

# Randomly initialize the factorization.
A <- matrix(runif(n*k),n,k)
B <- matrix(runif(m*k),k,m)

# RUN BASIC EM UPDATES
# --------------------
cat("Running the multiplcative (EM) updates.\n")
fit1 <- betanmf(X,A,B,200)

# RUN ACCELERATED EM UPDATES
# --------------------------
cat("Running the accelerated EM updates.\n")
fit2 <- daarbetanmf(X,A,B,200)

# RUN NNMF
# --------
cat("Running SQP updates.\n")
fit3 <- nnmf(X,k,init = list(W = A,H = B),method = "scd",loss = "mkl",
             rel.tol = 0,n.threads = 0,max.iter = 200,inner.max.iter = 4,
             trace = 1,verbose = 0)

# RUN ALT-SQP
# -----------
cat("Running extrapolated SQP updates.\n")
fit4 <- altsqp(X,list(L = A,F = t(B)),numiter = 200,
               control = list(extrapolate = 20),
               verbose = FALSE)

# PLOT IMPROVEMENT IN SOLUTION OVER TIME
# --------------------------------------
library(ggplot2)
library(cowplot)
f <- min(c(fit1$value,fit2$value,fit3$mkl,fit4$value)) - 0.01
pdat <-
  rbind(data.frame(iter = 1:200,dist = fit1$value - f,method = "EM"),
        data.frame(iter = 1:200,dist = fit2$value - f,method = "daarem"),
        data.frame(iter = 1:200,dist = fit3$mkl - f,method = "nnmf"),
        data.frame(iter = 1:200,dist = fit4$progress$objective - f,
                   method = "altsqp"))
p <- ggplot(pdat,aes(x = iter,y = dist,color = method)) +
  geom_line(size = 1) +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = c("darkblue","dodgerblue","magenta",
                                "darkorange")) +
  labs(x = "iteration",y = "distance from solution")
