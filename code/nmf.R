# TO DO: Explain here what this function does, and how to use it.
betanmf <- function (X, A, B, numiter = 100, e = 1e-15) {

  # Iterate the multiplicative (EM) updates.
  for (i in 1:numiter) {

    # Update the loadings ("activations").
    A <- scale.cols(A * ((X / (A %*% B)) %*% t(B)),1/rowSums(B))
    A <- pmax(A,e)

    # Update the factors ("basis vectors").
    B <- B * (t(A) %*% (X / (A %*% B))) / colSums(A)
    B <- pmax(B,e)

    # Record the algorithm's progress.
  }
}

# Scale each column A[,i] by b[i].
scale.cols <- function (A, b) 
  t(t(A) * b)
