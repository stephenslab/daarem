# A small script to illustrate application of the DAAREM method for
# computing maximum-likelihood estimates of mixture proportions in a
# mixture model.

# SET UP ENVIRONMENT
# ------------------
library(ggplot2)
library(cowplot)
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
out <- system.time(fit1 <- mixem(L,x0,numiter = 200))
f1  <- mixobjective(L,fit1$x)
cat(sprintf("Computation took %0.2f seconds.\n",out["elapsed"]))
cat(sprintf("Log-likelihood at EM estimate is %0.12f.\n",f1))

# RUN ACCELERATED EM
# ------------------
cat("Fitting mixture model with accelerated EM method (DAAREM).\n")
out <- system.time(fit2 <- mixdaarem(L,x0,numiter = 200))
f2  <- mixobjective(L,fit2$x)
cat(sprintf("Computation took %0.2f seconds.\n",out["elapsed"]))
cat(sprintf("Objective value at DAAREM estimate is %0.12f.\n",f2))

# PLOT IMPROVEMENT IN SOLUTION OVER TIME
# --------------------------------------
f    <- mixobjective(L,x)
pdat <-
  rbind(data.frame(iter = 1:200,dist = f - fit1$value,method = "EM"),
        data.frame(iter = 1:200,dist = f - fit2$value,method = "DAAREM"))
p <- ggplot(pdat,aes(x = iter,y = dist,col = method)) +
  geom_line(size = 1) +
  scale_y_continuous(trans = "log10",breaks = 10^seq(-4,4)) +
  scale_color_manual(values = c("darkorange","dodgerblue")) +
  labs(x = "iteration",y = "distance from solution")
print(p)
