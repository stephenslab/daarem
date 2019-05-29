# SET UP ENVIRONMENT
# ------------------
library(daarem)
source("misc.R")
source("mixem.R")

# LOAD DATA
# ---------
cat("Loading data.\n")
load("../data/mixdata.RData")
k <- ncol(L)

# RUN BASIC EM
# ------------
cat("Fitting mixture model with basic EM method.\n")
x0   <- rep(1/k,k)
out <- system.time(fit1 <- mixem(L,w,x0,numiter = 1000))
cat("Computation took %0.2f seconds.\n",out["elapsed"])

# RUN ACCELERATED EM
# ------------------
cat("Fitting mixture model with accelerated EM method.\n")
