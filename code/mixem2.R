bfgsmixem <- function (L, x0, numiter = 1000, e = 1e-15) {

  # Get the initial estimate of the solution.
  x <- log(x0)

  # Scale the correction factor (e) by the maximum value of each row
  # of the matrix (L). Therefore, we end up with a separate correction
  # factor for each row of L.
  e <- e * apply(L,1,max)

  # Initialize the BFGS approximation to the inverse of the Hessian
  # (which is denoted by "H" in Nocedal & Wright, 2nd ed).
  H <- diag(length(x))
  
  # This variable is used to keep track of the algorithm's progress;
  # it stores the value of the objective at each iteration.
  value <- rep(0,numiter)

  # Iterate the E and M steps.
  for (i in 1:numiter) {

    # Store the current iterate.
    x0 <- x

    # Obtain a descent direction.
    g <- x - log(mixem.update(L,softmax(x),e))

    # Compute the BFGS search direction.
    p <- solve(H,-g)
    
    # Run backtracking line search to find a suitable step size.
    a <- backtracking.line.search(x,g,p,
                                 function (x) -mixobjective(L,softmax(x),e),
                                 suffdecr = 0.01)
    
    # Move to the new iterate.
    x <- x + a*p
    print(a)

    # Update the BFGS approximation to the inverse Hessian.
    g0 <- g
    g  <- x - log(mixem.update(L,softmax(x),e))
    s  <- x - x0
    y  <- g - g0
    Hnew <- H + (y - H %*% s) %*% t(s)/dot(s,s)
    H <- H*0.9 + 0.1*Hnew
    
    # Record the algorithm's progress.
    value[i] <- mixobjective(L,softmax(x),e)
  }

  # Return the estimate of the solution ("x") and the value of the
  # objective at each iteration ("value").
  return(list(x = softmax(x),value = value))
}
