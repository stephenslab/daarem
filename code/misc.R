# Scale each column A[,i] by b[i].
scale.cols <- function (A, b)
  t(t(A) * b)

# Return the quadratic norm (2-norm) of vector x.
norm2 <- function (x)
  sqrt(sum(x*x))

# Compute the softmax of x.
softmax <- function (x) {
  y <- exp(c(0,x))
  return(y/sum(y))
}

# Compute the inverse softmax of y.
softmaxinverse <- function (y)
  log(y[-1]/y[1])
