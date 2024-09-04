# Install and load necessary package for truncated normal distribution
install.packages("truncnorm")
library(truncnorm)

# GCA Optimization Function
# This function performs optimization using the Ground and Air Forces Algorithm.
# Inputs:
# - f: Objective function to be minimized.
# - n: Total number of initial solutions (population size).
# - iter: Number of iterations.
# - bounds: A matrix with lower and upper bounds for each variable.
GCA <- function(f, n, iter, bounds) {
  # Extract lower and upper bounds for variables
  a <- bounds[, 1]
  b <- bounds[, 2]
  d <- nrow(bounds)  # Number of dimensions
  n <- n - d          # Adjust population size excluding air forces
  
  # Generate initial population within bounds
  X <- matrix(runif((n + d) * d, min = a, max = b), nrow = n + d, ncol = d)
  
  # Initialize personal best positions and their objective values
  pbest <- X
  pbestobj <- apply(pbest, 1, f)
  
  # Identify the global best position and value
  gbest_idx <- which.min(pbestobj)
  gbest <- pbest[gbest_idx, ]
  
  # Main optimization loop
  for (iteration in 1:iter) {
    # Update Ground Forces (GF)
    GFbest <- matrix(rep(gbest, n), nrow = n, ncol = d, byrow = TRUE)
    sigma_GF <- abs(GFbest - pbest[1:n, ])
    mu_GF <- (2 / 3) * GFbest + (1 / 3) * pbest[1:n, ]
    TN_GF <- matrix(rtruncnorm(n * d, a = a, b = b, mean = mu_GF, sd = sigma_GF), nrow = n, ncol = d)
    
    # Update Air Forces (AF)
    AFbest <- matrix(rep(gbest, d), nrow = d, ncol = d, byrow = FALSE)
    sigma_AF <- abs(AFbest - pbest[(n + 1):(n + d), ]) / 3
    mu_AF <- (2 / 3) * AFbest + (1 / 3) * pbest[(n + 1):(n + d), ]
    TN_AF <- matrix(rtruncnorm(d * d, a = a, b = b, mean = mu_AF, sd = sigma_AF), nrow = d, ncol = d)
    
    # Combine forces to create new population
    X <- rbind(TN_GF, TN_AF)
    
    # Evaluate new population
    obj <- apply(X, 1, f)
    
    # Update personal bests
    idx <- which(obj <= pbestobj)
    pbest[idx, ] <- X[idx, ]
    pbestobj[idx] <- obj[idx]
    
    # Update global best
    gbest_idx <- which.min(pbestobj)
    gbest <- pbest[gbest_idx, ]
  }
  
  # Print the best solution found
  print(gbest)
  z <- format(f(gbest), digits = 10)
  print(z)
}

# Example usage: Cantilever Beam Design Problem
# Define the objective function to be minimized
f <- function(x) {
  z <- 0.06224 * sum(x)
  t <- sum(c(61, 37, 19, 7, 1) / x^3)
  E <- (t - 1 > 0) * 100 * max(t - 1, 0) +
    (x - 100 > 0) * 100 * max(t - 100, 0) + (x - 0.01 < 0) * 100 * max(0.01 - x, 0)
  w <- z + E[1]
  return(w)
}

# Define bounds for the problem
bounds <- matrix(rep(c(0.01, 100), 5), nrow = 5, ncol = 2, byrow = TRUE)

# Run the GCA optimization
GCA(f, n = 100, iter = 500, bounds)
