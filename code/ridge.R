# TO DO: Explain here what this function does, and how to use it.
ridge <- function (X, y, s0, numiter = 100) {

  # Get the number of predictors.
  p <- ncol(X)

  # Compute some quantities used in the co-ordinate ascent updates below.
  xy <- drop(y %*% X)
  XX <- crossprod(X)
  d  <- diag(XX)
  
  # This variable is used to keep track of the algorithm's progress;
  # it stores the value of the objective at each iteration.
  value <- rep(0,numiter)
  
  # Iterate the co-ordinate ascent steps.
  for (iter in 1:numiter) {

    # Update the regression coefficients via co-ordinate ascent.
    for (i in 1:p)
      b[i] <- (xy[i] - sum(XX[i,-i] * b[-i]))/(d[i] + 1/s0)
        
    # Record the algorithm's progress.
    value[iter] <- ridge.objective(X,y,b,s0)
  }

  # Return the estimate of the regression coefficients ("b") and the
  # value of the objective at each iteration ("value").
  return(list(b = b,value = value))
}

# Compute the value of the log-posterior up to a constant of
# proportionality.
ridge.objective <- function (X, y, b, s0)
  -(norm2(y - X %*% b)^2 + norm2(b)^2/s0)
