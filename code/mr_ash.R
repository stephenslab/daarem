# TO DO: Explain here what this function does.
mr_ash <- function (X, y, b0, s, s0, w, numiter = 100) {

  # Get the initial estimates of the posterior mean regression coefficients.
  b <- b0
  
  # This variable is used to keep track of the algorithm's progress;
  # it stores the value of the objective (the variational lower bound,
  # or "ELBO") at each iteration.
  value <- rep(0,numiter)

  # Iterate the co-ordinate ascent steps.
  for (iter in 1:numiter) {

    # Run the co-ordinate ascent updates.
    b <- mr_ash_update(X,y,b,s,s0,w)

    # Record the algorithm's progress.
    value[iter] <- mr_ash_elbo(X,y,b,s,s0,w)
  }

  # Return the estimate of the regression coefficients ("b") and the
  # value of the objective at each iteration ("value").
  return(list(b = b,value = value))
}

# TO DO: Explain here what this function does.
daar_mr_ash <- function (X, y, b0, s, s0, w, numiter = 100, order = 10) {

  # Run DAAREM.
  out <- suppressWarnings(
    daarem(b0,daar_mr_ash_update,daar_mr_ash_objective,X,y,s,s0,w,
           control = list(maxiter = numiter,order = order,tol = 0,
                          mon.tol = 0.01,kappa = 20,alpha = 1.1)))

  # Return the estimate of the regression coefficients ("b") and the
  # value of the objective at each iteration ("value").
  return(list(b = out$par,value = out$objfn.track[-1]))
}

daar_mr_ash_update <- function (b, X, y, s, s0, w)
  mr_ash_update(X,y,b,s,s0,w)

daar_mr_ash_objective <- function (b, X, y, s, s0, w)
  mr_ash_elbo(X,y,b,s,s0,w)

# Compute the variational lower bound to the marginal log-likelihood.
mr_ash_elbo <- function (X, y, b, s, s0, w) {
  e   <- 1e-15
  k   <- length(w)
  d   <- colSums(X^2)
  out <- mr_ash_posterior(X,y,b,s,s0,w)
  a   <- out$a
  mu  <- out$mu
  v   <- out$v
  r   <- rowSums(a*mu)
  f   <- -norm2(y - X %*% r)^2/(2*s) - dot(d,betavarmix(a,mu,v))/(2*s)
  for (i in 1:k)
    f <- f + sum(a[,i]*log(w[i])) - sum(a[,i]*log(a[,i] + e)) +
             sum(a[,i])/2 + dot(a[,i],log(v[,i]/(s*s0[i])))/2 -
             dot(a[,i],v[,i] + mu[,i]^2)/(2*s*s0[i])
  return(f)
}

# Apply a co-ordinate ascent update once for each variable.
mr_ash_update <- function (X, y, b, s, s0, w) {

  # Get the number of variables (p) and the number of mixture
  # components (k).
  p <- length(b)
  k <- length(w)
  
  # Precompute some quantities used in the co-ordinate ascent updates
  # below.
  d  <- colSums(X^2)
  xy <- drop(y %*% X)

  # Compute the vector of posterior mean regression coefficients, r,
  # and the matrix-vector product X*r.
  out <- mr_ash_posterior(X,y,b,s,s0,w)
  r   <- rowSums(out$a * out$mu)
  xr  <- drop(X %*% r)

  # Repeat for each variable.
  for (i in 1:p) {
    r0 <- r[i]

    # Compute "b" using the ridge co-ordinate ascent update, in which
    # the prior on the coefficients is normal with zero mean and
    # variance 1.
    v    <- s/(d[i] + 1/s0)
    b[i] <- v[k]/s*(xy[i] + d[i]*r[i] - dot(xr,X[,i]))

    # Compute the new posterior mean regression coefficient.
    mu   <- v/v[k]*b[i]
    a    <- normalizelogweights(log(w) + (log(v/(s*s0)) + mu^2/v)/2)
    r[i] <- dot(a,mu)

    # Update the matrix-vector product X*r, where r is the vector of
    # posterior mean coefficients.
    xr <- xr + (r[i] - r0)*X[,i]
  }
  
  return(b)
}

# TO DO: Explain here what this function does.
mr_ash_posterior <- function (X, y, b, s, s0, w) {
  p  <- length(b)
  k  <- length(w)
  d  <- colSums(X^2)
  a  <- matrix(0,p,k)
  mu <- matrix(0,p,k)
  v  <- matrix(0,p,k)
  for (i in 1:p) {
    v[i,]  <- s/(d[i] + 1/s0)
    mu[i,] <- v[i,]/v[i,k]*b[i]
    a[i,]  <- normalizelogweights(log(w) + (log(v[i,]/(s*s0)) +
                                            mu[i,]^2/v[i,])/2)
  }
  return(list(a = a,mu = mu,v = v))
}

# Returns variances of variables drawn from mixtures of normals, in
# which vectors mu and s give the normal mean and variances, and
# vector p gives the mixture weights.
betavarmix <- function (p, mu, s)
  rowSums(p*(s + mu^2)) - rowSums(p*mu)^2
