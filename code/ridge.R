# Fit a ridge regression model by running a fixed number of
# co-ordinate ascent updates.
ridge <- function (X, y, b0, s0, numiter = 100) {

  # Get the initial estimate.
  b <- b0

  # This variable is used to keep track of the algorithm's progress;
  # it stores the value of the objective at each iteration.
  value <- rep(0,numiter)

  # Iterate the co-ordinate ascent steps.
  for (iter in 1:numiter) {

    # Run the co-ordinate ascent updates.
    b <- ridge.update(X,y,b,s0)
          
    # Record the algorithm's progress.
    value[iter] <- ridge.objective(X,y,b,s0)
  }
  
  # Return the estimate of the regression coefficients ("b") and the
  # value of the objective at each iteration ("value").
  return(list(b = b,value = value))
}

# Fit a ridge regression model by running a fixed number of
# co-ordinate ascent updates accelerated using the DAAREM method.
daarridge <- function (X, y, b0, s0, numiter = 100, order = 10) {
  nupd <- 0
  nobj <- 0
  
  # This implements the objfn argument for the daarem call below.
  daarridge.objective <- function (b, X, y, s0) {
    nobj <<- nobj + 1
    return(ridge.objective(X,y,b,s0))
  }
   # This implements the fixptfn argument for the daarem call below.
  daarridge.update <- function (b, X, y, x0) {
    nupd <<- nupd + 1
    return(ridge.update(X,y,b,s0))
  }
    
  # Run DAAREM.
  out <- suppressWarnings(
    daarem(b0,daarridge.update,daarridge.objective,X,y,s0,
           control = list(maxiter = numiter,order = order,tol = 0,
                          mon.tol = 0,kappa = 20,alpha = 1.2)))

  # Return the estimate of the regression coefficients ("b") and the
  # value of the objective at each iteration ("value").
  return(list(b     = out$par,
              value = out$objfn.track[-1],
              nupd  = nupd,
              nobj  = nobj))
}

# Compute the value of the log-posterior (i.e., the log-likelihood
# with the ridge "penalty term") up to a proportionality constant.
ridge.objective <- function (X, y, b, s0)
  -(norm2(y - X %*% b)^2 + norm2(b)^2/s0)

# Apply a co-ordinate ascent update once for each variable.
ridge.update <- function (X, y, b, s0) {
  p  <- length(b)
  d  <- colSums(X^2)
  xy <- drop(y %*% X)
  xb <- drop(X %*% b)
  for (i in 1:p) {
    b0   <- b[i]
    b[i] <- (xy[i] + d[i]*b[i] - dot(xb,X[,i]))/(d[i] + 1/s0)
    xb   <- xb + (b[i] - b0)*X[,i]
  }
  return(b)
}
