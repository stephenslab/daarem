# TO DO: Explain here what this function does, and how to use it.
betanmf <- function (X, A, B, numiter = 100, e = 1e-15) {

  # This variable is used to keep track of the algorithm's progress;
  # it stores the value of the cost function at each iteration.
  value <- rep(0,numiter)
  
  # Iterate the multiplicative (EM) updates.
  for (i in 1:numiter) {

    # Update the loadings ("activations").
    A <- scale.cols(A * ((X / (A %*% B)) %*% t(B)),1/rowSums(B))
    A <- pmax(A,e)

    # Update the factors ("basis vectors").
    B <- B * (t(A) %*% (X / (A %*% B))) / colSums(A)
    B <- pmax(B,e)

    # Record the algorithm's progress.
    AB       <- A %*% B
    value[i] <- sum(AB - X*log(AB + e))

  }

  # Return the estimated factors and loadings, and the value of the
  # cost function at each iteration ("value").
  return(list(A = A,B = B,value = value))
}

# Scale each column A[,i] by b[i].
scale.cols <- function (A, b) 
  t(t(A) * b)
