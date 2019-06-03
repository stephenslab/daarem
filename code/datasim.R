# Return a matrix of samples drawn from the multivariate normal
# distribution with zero mean and covariance S. The p x p covariance
# matrix S has entries S[i,j] = s^abs(i - j) for all i, j.
#
# Setting s = 0.5 reproduces the simulation of the predictors used in
# Examples 1 and 2 of Zou & Hastie (2005).
simulate_predictors_decaying_corr <- function (n, p, s = 0.5) {
  S <- s^abs(outer(1:p,1:p,"-"))
  return(mvrnorm(n,rep(0,p),S))
}
