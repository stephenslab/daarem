# TO DO: Explain here what this function does, and how to use it.
ridge <- function (X, y, b0, s0, numiter = 100) {

  # Get the initial estimate of the solution.
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

# TO DO: Explain here what this function does, and how to use it.
daarridge <- function (X, y, b0, s0, numiter = 100, order = 10) {

  # Get the initial estimate of the solution.
  b <- b0
  
  # Run DAAREM.
  out <- suppressWarnings(
    daarem(b,daarridge.update,daarridge.objective,X,y,s0,
           control = list(maxiter = numiter,order = order,tol = 0,
                          mon.tol = 0,kappa = 20,alpha = 1.2)))

  # Return the estimate of the regression coefficients ("b") and the
  # value of the objective at each iteration ("value").
  return(list(b = out$par,value = out$objfn.track[-1]))
}

# Compute the value of the log-posterior (i.e., the log-likelihood
# with the ridge "penalty term") up to a constant of proportionality.
ridge.objective <- function (X, y, b, s0)
  -(norm2(y - X %*% b)^2 + norm2(b/s0)^2)

# Run each of the co-ordinate ascent updates once for each co-ordinate.
ridge.update <- function (X, y, b, s0) {
  p <- length(b)
  xy <- drop(y %*% X)
  R  <- crossprod(X)
  for (i in 1:p)
    b[i] <- (xy[i] - sum(R[i,-i] * b[-i]))/(R[i,i] + 1/s0^2)
  return(b)
}

# This implements the objfn argument for the daarem call above.
daarridge.objective <- function (b, X, y, s0)
  ridge.objective(X,y,b,s0)

# This implements the fixptfn argument for the daarem call above.
daarridge.update <- function (b, X, y, x0)
  ridge.update(X,y,b,s0)
