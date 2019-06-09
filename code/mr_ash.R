# Fit a ridge regression model by running a fixed number of
# co-ordinate ascent updates.
mr_ash <- function (X, y, b0, s0, numiter = 100) {

  # Get the initial estimate.
  b <- b0
  
  # This variable is used to keep track of the algorithm's progress;
  # it stores the value of the objective at each iteration.
  value <- rep(0,numiter)
  
  # Iterate the co-ordinate ascent steps.
  for (iter in 1:numiter) {

    # Run the co-ordinate ascent updates.
    b <- ridge.update(X,y,b,s0)
          
    # Record the algorithm's progress.
    value[iter] <- ridge.objective(X,y,b,s0)
  }

  # Return the estimate of the regression coefficients ("b") and the
  # value of the objective at each iteration ("value").
  return(list(b = b,value = value))
}
