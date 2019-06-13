# Scale each column A[,i] by b[i].
scale.cols <- function (A, b)
  t(t(A) * b)

# Return the dot product of vectors x and y.
dot <- function (x,y)
  sum(x*y)

# Return the quadratic norm (2-norm) of vector x.
norm2 <- function (x)
  sqrt(dot(x,x))

# Compute the softmax of x.
softmax <- function (x) {
  y <- exp(x)
  return(y/sum(y))
}

# Compute the inverse softmax of y.
softmaxinverse <- function (y)
  log(y)

# Takes as input an array of unnormalized log-probabilities logw and
# returns normalized probabilities such that the sum is equal to 1.
normalizelogweights <- function (logw) {

  # Guard against underflow or overflow by adjusting the
  # log-probabilities so that the largest probability is 1.
  c <- max(logw)
  w <- exp(logw - c)

  # Normalize the probabilities.
  return(w/sum(w))
}
