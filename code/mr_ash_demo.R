## ----sim-settings--------------------------------------------------------
n  <- 200 
p  <- 500
na <- 50
se <- 4

# mr-ash model settings
s0 <- c(0.01,0.5,1)^2/se
w  <- c(0.9,0.05,0.05)

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
X <- simulate_predictors_decaying_corr(n,p,s = 0.8)
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

fit1 <- ridge(X,y,b0,s0[3],numiter = 200)
plot(log10(max(fit1$value) - fit1$value[1:100]),pch = 20)

# Fit mr-ash model using the basic co-ordinate ascent updates.
fit2 <- mr_ash(X,y,b0,se,s0,w,numiter = 200)
plot(log10(max(fit2$value) - fit2$value[1:190]),pch = 20)
