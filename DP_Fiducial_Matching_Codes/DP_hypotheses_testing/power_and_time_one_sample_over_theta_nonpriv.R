#' Compute Empirical Power and Timing for Binomial Hypothesis Testing
#'
#' This function simulates the empirical power of a hypothesis test for binomial proportions
#' under the null (`theta0`) and alternative (`theta_values`) hypotheses. It also records
#' computation time for each simulation.
#'
#' @param theta0 Numeric: Null hypothesis probability (must be in [0, 1]).
#' @param theta_values Numeric vector: Alternative hypothesis probabilities to test (must be in [0, 1]).
#' @param n Integer: Sample size (must be > 0).
#' @param B Integer: Number of Monte Carlo replications (must be > 0).
#' @param alpha Numeric: Significance level (default = 0.05, must be in (0, 1)).
#'
#' @return A list containing:
#'   - `results`: A `data.frame` with columns `theta`, `power_nonpriv`, `mean_time_nonpriv`, `median_time_nonpriv`.
#'   - `p_values`: A matrix of p-values (B rows × `length(theta_values)` columns).
#' @examples
#' theta0 <- 0.2
#' theta_values <- seq(0.2, 0.9, by = 0.1)
#' n <- 100
#' B <- 100
#' results <- compute_power(theta0, theta_values, n, B)
compute_power <- function(theta0, theta_values, n, B, alpha = 0.05) {
  # Input validation
  stopifnot(
    theta0 >= 0 && theta0 <= 1,
    all(theta_values >= 0 & theta_values <= 1),
    n > 0,
    B > 0,
    alpha > 0 && alpha < 1
  )
  
  # Initialize storage for p-values and compute times
  p_values <- matrix(NA, nrow = B, ncol = length(theta_values))
  compute_times <- matrix(NA, nrow = B, ncol = length(theta_values))
  
  # Precompute possible binomial values (0 to n)
  values <- seq(0, n)
  pdf_H0 <- dbinom(values, size = n, prob = theta0)  # PDF under H0
  
  # Loop over theta values
  for (j in seq_along(theta_values)) {
    theta <- theta_values[j]
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    
    # Monte Carlo simulation
    for (i in 1:B) {
      set.seed(i)  # Ensure reproducibility
      start_time <- Sys.time()
      
      # Simulate data under alternative hypothesis (theta)
      X <- rbinom(n = 1, size = n, prob = theta)
      
      # Add uniform noise for randomized p-value
      Z <- X + runif(n = 1, min = -0.5, max = 0.5)
      
      # Compute p-value using CDF of uniform distribution
      cdf <- punif(values - Z, min = -0.5, max = 0.5)
      p_values[i, j] <- sum(cdf * pdf_H0)
      
      # Record computation time (milliseconds)
      compute_times[i, j] <- as.numeric(Sys.time() - start_time) * 1000
      setTxtProgressBar(pb, i)
    }
    close(pb)
    message(sprintf("Completed theta = %.2f", theta))
  }
  
  # Summarize results
  results <- data.frame(
    theta = theta_values,
    power_nonpriv = colMeans(p_values < alpha),
    mean_time_nonpriv = colMeans(compute_times),
    median_time_nonpriv = apply(compute_times, 2, median)
  )
  
  return(results)
}

# Simulation parameters
theta0 <- 0.2  # Null hypothesis
theta_values <- seq(0.2, 0.9, by = 0.1)  # Alternatives to test
n <- 100       # Sample size
B <- 10^4       # Number of Monte Carlo replications
alpha <- 0.05  # Significance level

# Run simulation
power_results <- compute_power(theta0, theta_values, n, B, alpha)
power_results 
