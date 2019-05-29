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
n <- nrow(L)
m <- ncol(L)
cat(sprintf("Loaded %d x %d data matrix.\n",n,m))

# Set initial estimate of mixture proportions.
x0 <- rep(1/m,m)

# RUN BASIC EM
# ------------
cat("Fitting mixture model with basic EM method.\n")
out <- system.time(fit1 <- mixem(L,x0,numiter = 1000))
f1  <- mixobjective(L,fit1$x)
cat(sprintf("Computation took %0.2f seconds.\n",out["elapsed"]))
cat(sprintf("Log-likelihood at EM estimate is %0.12f.\n",f1))

# RUN ACCELERATED EM
# ------------------
cat("Fitting mixture model with accelerated EM method (DAAREM).\n")
out <- system.time(fit2 <- mixdaarem(L,x0,numiter = 1000))
f2  <- mixobjective(L,fit2$x)
cat(sprintf("Computation took %0.2f seconds.\n",out["elapsed"]))
cat(sprintf("Objective value at DAAREM estimate is %0.12f.\n",f2))

# PLOT IMPROVEMENT IN SOLUTION OVER TIME
# --------------------------------------
# TO DO.
