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

  # This variable is used to keep track of the algorithm's progress;
  # it stores the value of the objective at each iteration.
  value <- rep(0,numiter)

  # Iterate the E and M steps.
  for (i in 1:numiter) {

    # Update the solution.
    x <- mixem.update(L,x,e)
    
    # Record the algorithm's progress.
    value[i] <- mixobjective(L,x,e)
  }

  # Return the estimate of the solution ("x") and the value of the
  # objective at each iteration ("value").
  return(list(x = x,value = value))
}

# Compute maximum-likelihood estimates of the mixture proportions in a
# mixture model by running DAAREM, an accelerated variant of the EM
# algorithm.
mixdaarem <- function (L, x0, numiter = 1000, order = 10, e = 1e-15) {

  # Get the initial estimate of the solution.
  x <- x0
    
  # Scale the correction factor (e) by the maximum value of each row
  # of the matrix (L). Therefore, we end up with a separate correction
  # factor for each row of L.
  e <- e * apply(L,1,max)

  # Run DAAREM.
  out <- suppressWarnings(
    daarem(x,mixdaarem.update,mixdaarem.objective,L,e,
           control = list(maxiter = numiter,order = order,tol = 0,
                          mon.tol = 0.05,kappa = 20,alpha = 1.2)))

  # Return the estimate of the solution and the value of the objective
  # at each iteration.
  return(list(x = project.iterate(out$par),value = out$objfn.track[-1]))
}

# Project the iterate so that it lies on the simplex.
project.iterate <- function (x) {
  x <- pmax(x,0)
  return(x/sum(x))
}

# This implements the fixptfn argument for the daarem call above.
mixdaarem.update <- function (x, L, e)
  mixem.update(L,project.iterate(x),e)

# This implements the objfn argument for the daarem call above.
mixdaarem.objective <- function (x, L, e)
  mixobjective(L,project.iterate(x),e)

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
