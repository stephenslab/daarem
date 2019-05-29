# Compute maximum-likelihood estimates of the mixture proportions in a
# mixture model by iterating the EM updates for a fixed number of
# iterations.
mixem <- function (L, x0, numiter = 1000, e = 1e-15) {

  # Get the initial estimate of the solution.
  x <- x0

  # Scale the correction factor (e) by the maximum value of each row
  # of the matrix (L). Therefore, we end up with a separate correction
  # factor for each row of L.
  e <- e * apply(L,1,max)

  # This data frame is used to store the EM algorithm's progress at
  # each iteration. The three columns of the data frame store: (1) the
  # iteration number, (2) the objective value at each iteation, and (3)
  # the largest change in the solution at each iteration.
  progress <- data.frame(iter = 1:numiter,obj = 0,maxd = 0)

  # Iterate the E and M steps.
  for (i in 1:numiter) {

    # Store the current estimate of the solution.
    x0 <- x

    # Update the solution.
    x <- mixem.update(L,x,e)
    
    # Record the algorithm's progress.
    progress[i,"obj"]  <- mixobjective(L,x,e)
    progress[i,"maxd"] <- max(abs(x - x0))
  }

  # Return (1) the estimate of the solution, (2) the value of the
  # objective at this estimate, and (3) a record of the progress made
  # at each EM iteration.
  f <- mixobjective(L,x,e)
  return(list(x = x,value = f,progress = progress))
}

# TO DO: Explain here what this function does, and how to use it.
mixdaarem <- function (L, x0, numiter = 1000, e = 1e-15) {

  # Get the initial estimate of the solution.
  x <- log(x0 + e)

  # Run DAAREM.
  out <- suppressWarnings(daarem(x,mixdaarem.update,mixdaarem.objective,L,e,
                                 control = list(maxiter = numiter,tol = 0)))

  return(out)
}

# This implements the fixptfn argument for the daarem call above.
mixdaarem.update <- function (x, L, e) {
  x <- softmax(x)
  x <- mixem.update(L,x,e)
  return(log(x))
}

# This implements the objfn argument for the daarem call above.
mixdaarem.objective <- function (x, L, e) {
  return(mixobjective(L,softmax(x),e))
}

# Perform a single expectation maximization (EM) update.
mixem.update <- function (L, x, e) {

  # E STEP
  # ------
  # Compute the n x m matrix of posterior mixture assignment
  # probabilities (L is an n x m matrix).
  P <- scale.cols(L,x) + e
  P <- P / rowSums(P)

  # M STEP
  # ------
  # Update the mixture weights.
  return(colMeans(P))
}

# Compute the value of the log-likelihood at x; e is a vector in which
# the entries can be set to small, positive numbers, or to zero.
mixobjective <- function (L, x, e = 1e-15) {
 y <- drop(L %*% x) + e
 if (all(y > 0))
   return(sum(log(y)))
 else
   return(-Inf)
}
