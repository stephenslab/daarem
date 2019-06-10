## ----sim-settings--------------------------------------------------------
n  <- 200 
p  <- 500
na <- 10
se <- 4

# mr-ash model settings
s0 <- c(0.1,1,10)^2/se
w  <- c(0.5,0.25,0.25)

## ----load-pkgs, warning=FALSE, message=FALSE-----------------------------
library(MASS)
library(daarem)
library(ggplot2)
library(cowplot)
source("../code/misc.R")
source("../code/datasim.R")
source("../code/ridge.R")
source("../code/mr_ash.R")

set.seed(1)

## ----sim-x---------------------------------------------------------------
X <- simulate_predictors_decaying_corr(n,p,s = 0.5)
X <- scale(X,center = TRUE,scale = FALSE)

## ----sim-beta------------------------------------------------------------
i    <- sample(p,na)
b    <- rep(0,p)
b[i] <- rnorm(na)
 
## ----sim-y---------------------------------------------------------------
y <- drop(X %*% b + sqrt(se)*rnorm(n))
y <- y - mean(y)

## ----init-estimate-------------------------------------------------------
b0 <- rep(0,p)

## ----ridge-solution------------------------------------------------------
bhat <- drop(solve(t(X) %*% X + diag(rep(se^2,p)),t(X) %*% y))

# Fit mr-ash model using the basic co-ordinate ascent updates.
fit2 <- mr_ash(X,y,b0,se,s0,w,numiter = 200)
fit3 <- daar_mr_ash(X,y,b0,se,s0,w,numiter = 200)
fbest <- max(c(fit2$value,fit3$value))
plot(log10(fbest - fit2$value),pch = 20)
plot(log10(fbest - fit3$value),pch = 20)
