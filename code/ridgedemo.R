# TO DO: Explain here what this script does, and how to use it.

# SCRIPT PARAMETERS
# -----------------
# TO DO: Explain here what these parameters are for.
n  <- 200 
p  <- 500
na <- 10
se <- 2
s0 <- 1/se

# SET UP ENVIRONMENT
# ------------------
library(MASS)
library(daarem)
library(ggplot2)
library(cowplot)
source("misc.R")
source("datasim.R")
source("ridge.R")

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# SIMULATE DATA
# -------------
# Simulate the predictors.
cat("Simulating data.\n")
X <- simulate_predictors_decaying_corr(n,p,s = 0.5)
X <- scale(X,center = TRUE,scale = FALSE)

# Generate additive effects for the markers so that exactly na of them have
# a nonzero effect on the trait.
i    <- sample(p,na)
b    <- rep(0,p)
b[i] <- rnorm(na)

# Simulate the continuous outcomes.
y <- drop(X %*% b + se*rnorm(n))
y <- y - mean(y)

# Set initial estimate of mixture proportions.
b0 <-  

# FIT RIDGE REGRESSION MODEL
# --------------------------
# TO DO: Explain here what these lines of code do.
cat("Fitting ridge regression via basic co-ordinate ascent updates.\n")
out <- system.time(fit1 <- ridge(X,y,s0))
f1  <- ridge.objective(X,y,fit1$b,s0)
cat(sprintf("Computation took %0.2f seconds.\n",out["elapsed"]))
cat(sprintf("Objective value at solution is %0.12f.\n",f1))
