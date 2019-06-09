## ----knitr, echo=FALSE---------------------------------------------------
knitr::opts_chunk$set(comment = "#",results = "hold",collapse = TRUE,
                      fig.align = "center")

## ----sim-settings--------------------------------------------------------
n  <- 200 
p  <- 500
na <- 10
se <- 2

## ----model-settings------------------------------------------------------
s0 <- 1/se

## ----load-pkgs, warning=FALSE, message=FALSE-----------------------------
library(MASS)
library(daarem)
library(ggplot2)
library(cowplot)
source("../code/misc.R")
source("../code/datasim.R")
source("../code/ridge.R")

## set.seed(1)

## ----sim-x---------------------------------------------------------------
X <- simulate_predictors_decaying_corr(n,p,s = 0.5)
X <- scale(X,center = TRUE,scale = FALSE)

## ----sim-beta------------------------------------------------------------
i    <- sample(p,na)
b    <- rep(0,p)
b[i] <- rnorm(na)

## ----sim-y---------------------------------------------------------------
y <- drop(X %*% b + se*rnorm(n))
y <- y - mean(y)

## ----init-estimate-------------------------------------------------------
b0 <- rep(0,p)

## ----fit-ridge-----------------------------------------------------------
out <- system.time(fit1 <- ridge(X,y,b0,s0,numiter = 100))
f1  <- ridge.objective(X,y,fit1$b,s0)
cat(sprintf("Computation took %0.2f seconds.\n",out["elapsed"]))
cat(sprintf("Objective value at solution is %0.12f.\n",f1))

## ----fit-ridge-2---------------------------------------------------------
out <- system.time(fit2 <- daarridge(X,y,b0,s0,numiter = 100))
f2  <- ridge.objective(X,y,fit2$b,s0)
cat(sprintf("Computation took %0.2f seconds.\n",out["elapsed"]))
cat(sprintf("Objective value at solution is %0.12f.\n",f2))

## ----ridge-solution------------------------------------------------------
bhat <- drop(solve(t(X) %*% X + diag(rep(1/s0^2,p)),t(X) %*% y))
f    <- ridge.objective(X,y,bhat,s0)

## ----plot-iter-vs-objective, fig.height=4, fig.width=6-------------------
pdat <-
  rbind(data.frame(iter = 1:100,dist = f - fit1$value,method = "basic"),
        data.frame(iter = 1:100,dist = f - fit2$value,method = "accelerated"))
p <- ggplot(pdat,aes(x = iter,y = dist,col = method)) +
  geom_line(size = 1) +
  scale_y_continuous(trans = "log10",breaks = 10^seq(-8,4)) +
  scale_color_manual(values = c("darkorange","dodgerblue")) +
  labs(x = "iteration",y = "distance from solution")
print(p)

