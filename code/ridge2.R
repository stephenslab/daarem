# TO DO: Explain here what this function does, and how to use it.
bfgsridge <- function (X, y, b0, s0, numiter = 100) {

  # Get the initial estimate.
  b <- b0

  # Get the number of coefficients to estimate.
  n <- length(b)
  
  # Initialize the BFGS approximation to the inverse of the Hessian
  # (which is denoted by "H" in Nocedal & Wright, 2nd ed).
  H <- diag(n)
  
  # This variable is used to keep track of the algorithm's progress;
  # it stores the value of the objective at each iteration.
  value <- rep(0,numiter)

  # Iterate the coordinate-wise updates.
  for (iter in 1:numiter) {

    # Store the current iterate.
    b0 <- b

    # Run the coordinate-wise updates to obtain a descent direction.
    g <- b - ridge.update(X,y,b,s0)

    # Compute the BFGS search direction.
    p <- solve(H,-g)

    # Run backtracking line search to find a suitable step size.
    a <- backtracking.line.search(b,g,p,
                                  function (b) -ridge.objective(X,y,b,s0))
    
    # Move to the new iterate.
    b <- b + a*p
    
    # Update the BFGS approximation to the inverse Hessian.
    g0 <- g
    g  <- b - ridge.update(X,y,b,s0)
    s  <- b - b0
    u  <- g - g0
    H  <- H + (u - H %*% s) %*% t(s)/dot(s,s)
    
    # g0 <- b0/s0 - drop(drop(y - X %*% b0) %*% X)
    # g  <- b/s0 - drop(drop(y - X %*% b) %*% X)
    # H  <- bfgs.update(H,b - b0,g - g0)
    # print(dot(b - b0,g - g0))
    
    # Record the algorithm's progress.
    value[iter] <- ridge.objective(X,y,b,s0)
  }

  # Return the estimate of the regression coefficients ("b") and the
  # value of the objective at each iteration ("value").
  return(list(b = b,value = value))
}
