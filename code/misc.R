# Scale each column A[,i] by b[i].
scale.cols <- function (A, b)
  t(t(A) * b)

# Return the softmax of x.
softmax <- function (x)
  exp(x)/sum(exp(x))
