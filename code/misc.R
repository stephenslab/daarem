# Scale each column A[,i] by b[i].
scale.cols <- function (A, b)
  t(t(A) * b)

# Return the quadratic norm (2-norm) of vector x.
norm2 <- function (x)
  sqrt(sum(x*x))
