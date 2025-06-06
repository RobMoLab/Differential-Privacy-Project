#' Compute Empirical Power for Exact Binomial Test
#'
#' Simulates the power of an exact binomial test comparing a null hypothesis (theta0)
#' against alternative values (theta_values) using Monte Carlo methods.
#'
#' @param theta0 Numeric: Null hypothesis probability (must be in [0, 1]).
#' @param theta_values Numeric vector: Alternative probabilities to test (must be in [0, 1]).
#' @param n Integer: Sample size (must be > 0).
#' @param B Integer: Number of Monte Carlo replications (must be > 0).
#' @param alpha Numeric: Significance level (default = 0.05, must be in (0, 1)).
#'
#' @return A data.frame with columns:
#'   - theta: Alternative hypothesis values
#'   - power_exact: Empirical power of exact binomial test
#'   - mean_time_exact: Mean computation time in milliseconds
#' @examples
#' theta0 <- 0.2
#' theta_values <- seq(0.2, 0.9, by = 0.1)
#' n <- 100
#' B <- 1000
#' results <- compute_power_normal_approx(theta0, theta_values, n, B)
compute_power_normal_approx <- function(theta0, theta_values, n, B, alpha = 0.05) {
  # Input validation
  stopifnot(
    theta0 >= 0 && theta0 <= 1,
    all(theta_values >= 0 & theta_values <= 1),
    n > 0,
    B > 0,
    alpha > 0 && alpha < 1
  )
  
  # Initialize results dataframe
  results <- data.frame(
    theta = theta_values,
    power_exact = numeric(length(theta_values)),
    mean_time_exact = numeric(length(theta_values))
  )
  
  # Main simulation loop
  for (j in seq_along(theta_values)) {
    p_values_exact <- numeric(B)
    times_exact <- numeric(B)
    
    # Progress bar setup
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    message(sprintf("\nRunning simulations for theta = %.1f...", theta_values[j]))
    
    for (i in 1:B) {
      set.seed(i)  # Ensures reproducibility
      
      # --- Exact Binomial Test ---
      start_time <- Sys.time()
      x <- rbinom(1, size = n, prob = theta_values[j])
      p_values_exact[i] <- binom.test(x, n, p = theta0, 
                                      alternative = "greater")$p.value
      times_exact[i] <- as.numeric(Sys.time() - start_time) * 1000  # convert time to ms
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Store results
    results$power_exact[j] <- mean(p_values_exact < alpha)
    results$mean_time_exact[j] <- mean(times_exact)
    
    message(sprintf("Theta = %.1f: Power = %.3f, Mean Time = %.2f ms",
                    theta_values[j], results$power_exact[j], results$mean_time_exact[j]))
  }
  
  return(results)
}

# Simulation parameters
theta0 <- 0.2  # Null hypothesis probability
theta_values <- seq(0.2, 0.9, by = 0.1)  # Test values
n <- 100        # Sample size
B <- 10^4      # Replications (large B reduces Monte Carlo error)
alpha <- 0.05  # Significance level

# Run simulation with reproducibility
power_results <- compute_power_normal_approx(theta0, theta_values, n, B, alpha)

