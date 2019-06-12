# TO DO: Explain here what this function does, and how to use it.
#
# This uses equation 6.17 from Nocedal & Wright, 2nd ed.
bfgs.update <- function (H, s, y) {
  n <- length(s)
  r <- dot(s,y)
  if (r > 1e-6) {
    U <- diag(n) - s %*% t(y)/r
    return(U %*% H %*% t(U) + s %*% t(s)/r)
  } else
    return(H)
}

# TO DO: Explain here what this function does, and how to use it.
#
# NOTE: In the longer run, I will need a line search method that
# better ensures the curvature condition; see the discussion of
# backtrackinf line search in Nocedal & Wright, 2nd ed.
backtracking.line.search <- function (x, g, p, objfunc, suffdecr = 0.01,
                                      stepsizereduce = 0.75,
                                      minstepsize = 1e-8) {
  a <- 1
  f <- objfunc(x)
  while (TRUE) {
    y    <- x + a*p
    fnew <- objfunc(y)

    # Check whether the new candidate solution satisfies the
    # sufficient decrease condition.
    if (fnew < f + suffdecr*a*dot(p,g))
      break

    # If we cannot decrease the step size further, terminate the
    # backtracking line search, and set the step size to be the
    # minimum step size.
    else if (a*stepsizereduce < minstepsize) {
      y <- x + minstepsize*p
      break
    }
    
    # The new candidate does not satisfy the sufficient decrease
    # condition, so we need to try again with a smaller step size.
    a <- a * stepsizereduce
  }

  return(a)
}
