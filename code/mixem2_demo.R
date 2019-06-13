library(ggplot2)
library(cowplot)
library(daarem)
source("misc.R")
source("bfgs.R")
source("mixem.R")
source("mixem2.R")

## ----load-data-----------------------------------------------------------
load("../data/mixdata.RData")
n <- nrow(L)
m <- ncol(L)
cat(sprintf("Loaded %d x %d data matrix.\n",n,m))

## ----init-estimate-------------------------------------------------------
x0 <- rep(1/m,m)

## ----fit-em--------------------------------------------------------------
out <- system.time(fit1 <- mixem(L,x0,numiter = 200))
f1  <- mixobjective(L,fit1$x)
cat(sprintf("Computation took %0.2f seconds.\n",out["elapsed"]))
cat(sprintf("Log-likelihood at EM estimate is %0.12f.\n",f1))

## ----fit-daarem----------------------------------------------------------
out <- system.time(fit2 <- mixdaarem(L,x0,numiter = 200))
f2  <- mixobjective(L,fit2$x)
cat(sprintf("Computation took %0.2f seconds.\n",out["elapsed"]))
cat(sprintf("Objective value at DAAREM estimate is %0.12f.\n",f2))

## ----fit-bfgs------------------------------------------------------------
fit3 <- bfgsmixem(L,x0,numiter = 200)
f3  <- mixobjective(L,fit3$x)
cat(sprintf("Objective value at quasi-Newton estimate is %0.12f.\n",f2))

## ----plot-iter-vs-objective, fig.height=4, fig.width=6-------------------
f    <- mixobjective(L,x)
pdat <-
  rbind(data.frame(iter = 1:200,dist = f - fit1$value,method = "EM"),
        data.frame(iter = 1:200,dist = f - fit2$value,method = "DAAREM"),
        data.frame(iter = 1:200,dist = f - fit3$value,method = "Broyden"))
p <- ggplot(pdat,aes(x = iter,y = dist,col = method)) +
  geom_line(size = 1) +
  scale_y_continuous(trans = "log10",breaks = 10^seq(-4,4)) +
  scale_color_manual(values = c("darkorange","dodgerblue","magenta")) +
  labs(x = "iteration",y = "distance from solution")
print(p)
