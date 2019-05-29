# SET UP ENVIRONMENT
# ------------------
library(daarem)

# LOAD DATA
# ---------
cat("Loading data.\n")

# RUN BASIC EM
# ------------
cat("Fitting mixture model with basic EM method.\n")
fit1 <- mixem()

# RUN ACCELERATED EM
# ------------------
cat("Fitting mixture model with accelerated EM method.\n")
