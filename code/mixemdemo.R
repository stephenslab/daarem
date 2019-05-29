# TO DO: Briefly explain here what this script does, and how to use it.

# SET UP ENVIRONMENT
# ------------------
library(daarem)
source("misc.R")
source("mixem.R")

# LOAD DATA
# ---------
cat("Loading data.\n")
load("../data/mixdata.RData")

# Set initial estimate of mixture proportions.
k  <- ncol(L)
x0 <- rep(1/k,k)

# RUN BASIC EM
# ------------
cat("Fitting mixture model with basic EM method.\n")
out <- system.time(fit1 <- mixem(L,x0,numiter = 1000))
cat(sprintf("Computation took %0.2f seconds.\n",out["elapsed"]))
cat(sprintf("Log-likelihood at EM estimate is %0.12f.\n",fit1$value))

stop()

# RUN ACCELERATED EM
# ------------------
cat("Fitting mixture model with accelerated EM method.\n")
# TO DO.
cat(sprintf("Computation took %0.2f seconds.\n",out["elapsed"]))
cat(sprintf("Objective value at DAAREM estimate is %0.12f.\n",fit2$value))

# PLOT IMPROVEMENT IN SOLUTION OVER TIME
# --------------------------------------
# TO DO.
