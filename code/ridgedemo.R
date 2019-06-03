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
b0 <- rep(0,p)

# RUN BASIC CO-ORDINATE ASCENT UPDATES
# ------------------------------------
# TO DO: Explain here what these lines of code do.
cat("Fitting ridge regression with basic co-ordinate ascent updates.\n")
out <- system.time(fit1 <- ridge(X,y,b0,s0,numiter = 150))
f1  <- ridge.objective(X,y,fit1$b,s0)
cat(sprintf("Computation took %0.2f seconds.\n",out["elapsed"]))
cat(sprintf("Objective value at solution is %0.12f.\n",f1))

# RUN ACCELERATED CO-ORDINATE ASCENT UPDATES
cat("Fitting ridge regression with accelerated co-ordinate ascent updates.\n")
out <- system.time(fit2 <- daarridge(X,y,b0,s0,numiter = 150))
f2  <- ridge.objective(X,y,fit2$b,s0)
cat(sprintf("Computation took %0.2f seconds.\n",out["elapsed"]))
cat(sprintf("Objective value at solution is %0.12f.\n",f2))

# PLOT IMPROVEMENT IN SOLUTION OVER TIME
# --------------------------------------
bhat <- drop(solve(t(X) %*% X + diag(p)/s0,t(X) %*% y))
f    <- ridge.objective(X,y,bhat,s0)
pdat <-
  rbind(data.frame(iter = 1:150,dist = f - fit1$value,method = "basic"),
        data.frame(iter = 1:150,dist = f - fit2$value,method = "accelerated"))
p <- ggplot(pdat,aes(x = iter,y = dist,col = method)) +
  geom_line(size = 1) +
  scale_y_continuous(trans = "log10",breaks = 10^seq(-8,4)) +
  scale_color_manual(values = c("darkorange","dodgerblue")) +
  labs(x = "iteration",y = "distance from solution")
print(p)
