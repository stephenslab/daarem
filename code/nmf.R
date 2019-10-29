# TO DO: Explain here what this function does, and how to use it.
betanmf <- function (X, A, B, numiter = 100, e = 1e-15) {

  # This variable is used to keep track of the algorithm's progress;
  # it stores the value of the cost function at each iteration.
  value <- rep(0,numiter)
  
  # Iterate the multiplicative (EM) updates.
  cat("iter objective (cost fn)\n")
  for (i in 1:numiter) {
    out      <- betanmf.update(X,A,B,e)
    A        <- out$A
    B        <- out$B
    value[i] <- cost(X,A,B,e)
    cat(sprintf("%4d %+0.12e\n",i,value[i]))
  }

  # Return the estimated factors and loadings, and the value of the
  # cost function at each iteration ("value").
  return(list(A = A,B = B,value = value))
}

# TO DO: Explain here what this function does, and how to use it.
daarbetanmf <- function (X, A, B, numiter = 100,order = 10,e = 1e-15) {

  # Run DAAREM.
  out <- suppressWarnings(
    daarem(x0,daarbetanmf.update,daarbetanmf.objective,X,e,
           control = list(maxiter = numiter,order = order,tol = 0,
                          mon.tol = 0.05,kappa = 20,alpha = 1.2)))

  # Return the estimate of the solution and the value of the objective
  # at each iteration.
  # TO DO.
}

# This implements the fixptfn argument for the daarem call above.
daarbetanmf.update <- function (vars, X, e) {
  # TO DO.
}

# This implements the objfn argument for the daarem call above.
daarbetanmf.objective <- function (vars, X, e) {
  # TO DO.
}

# Perform a single multiplicative (EM) update.
betanmf.update <- function (X, A, B, e) {

  # Update the loadings ("activations").
  A <- scale.cols(A * ((X / (A %*% B)) %*% t(B)),1/rowSums(B))
  A <- pmax(A,e)

  # Update the factors ("basis vectors").
  B <- B * (t(A) %*% (X / (A %*% B))) / colSums(A)
  B <- pmax(B,e)

  return(list(A = A,B = B))
}

# Compute the value of the cost function for non-negative matrix
# factorization, in which matrix X is approximated by matrix AB = A*B.
# This is equivalent to the negative Poisson log-likelihood after
# removing terms that do not depend on A or B.
cost <- function (X, A, B, e) {
  AB <- A %*% B
  return(sum(AB - X*log(AB + e)))
}

# Scale each column A[,i] by b[i].
scale.cols <- function (A, b) 
  t(t(A) * b)

