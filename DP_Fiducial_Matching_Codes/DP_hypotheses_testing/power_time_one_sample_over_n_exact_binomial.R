#' Compute Power of Exact Binomial Test
#' 
#' Simulates the power of an exact binomial test for H₀: θ ≤ θ₀ vs Hₐ: θ > θ₀
#' across different sample sizes.
#' 
#' @param theta0 The threshold value for the null hypothesis (θ₀)
#' @param theta_true The true proportion under the alternative hypothesis (θ)
#' @param n_vec Vector of sample sizes to evaluate
#' @param B Number of Monte Carlo replications per sample size
#' @param alpha Significance level
#' 
#' @return A data frame containing:
#'   - n: Sample sizes
#'   - power: Estimated power at each sample size
#'   - mean_time_ms: Average computation time in milliseconds
#' 
#' @details
#' For each sample size in n_vec:
#' 1. Generates B binomial samples from Binom(n, θ)
#' 2. Performs exact binomial test (one-sided greater than)
#' 3. Computes proportion of p-values < α (power)
#' 4. Tracks computation times
#' 
#' Includes progress bar and interim results messages.
compute_exact_power <- function(theta0, theta_true, n_vec, B, alpha) {
  # Initialize results data frame
  results <- data.frame(
    n = n_vec,
    power = numeric(length(n_vec)),
    mean_time_ms = numeric(length(n_vec))
  )
  
  # Loop through sample sizes
  for (j in seq_along(n_vec)) {
    p_values <- numeric(B)  # Store p-values
    compute_times <- numeric(B)  # Store computation times
    
    # Initialize progress bar
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    
    # Monte Carlo replications
    for (i in 1:B) {
      set.seed(i)  # Reproducibility
      start_time <- Sys.time()
      
      # Generate data from alternative distribution
      X <- rbinom(1, size = n_vec[j], prob = theta_true)
      
      # Perform exact binomial test (one-sided greater than)
      test_result <- binom.test(X, n_vec[j], p = theta0, alternative = "greater")
      p_values[i] <- test_result$p.value
      
      # Record computation time (milliseconds)
      compute_times[i] <- as.numeric(Sys.time() - start_time) * 1000
      
      # Update progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)  # Close progress bar
    
    # Calculate power and average time
    results$power[j] <- mean(p_values < alpha)
    results$mean_time_ms[j] <- mean(compute_times)
    
    # Print interim results
    message(sprintf("n = %d: Power = %.3f | Time = %.2f ms", 
                    n_vec[j], results$power[j], results$mean_time_ms[j]))
  }
  
  return(results)
}


# Simulation parameters
n_vec <- c(16, 30, 50, 100, 150, 200, 350, 400, 500)  # Sample sizes to evaluate
theta0 <- 0.9      # Null hypothesis threshold
theta_true <- 0.95 # True proportion under alternative
B <- 10^4          # Number of Monte Carlo replications
alpha <- 0.05      # Significance level

# Run power simulation
results <- compute_exact_power(theta0, theta_true, n_vec, B, alpha)

# View results
print(results)

# Optional: Save results
# write.csv(results, "binomial_test_power_results.csv", row.names = FALSE)
