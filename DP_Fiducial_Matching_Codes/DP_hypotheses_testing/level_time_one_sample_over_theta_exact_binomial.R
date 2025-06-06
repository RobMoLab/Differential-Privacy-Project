#' Compute Exact Binomial Test Properties
#' 
#' Evaluates Type I error rates and computation times for exact binomial tests
#' across different true proportions (θ₀).
#' 
#' @param theta0 Vector of true proportions to evaluate
#' @param n Sample size (fixed across simulations)
#' @param B Number of Monte Carlo replications per θ₀ value
#' @param alpha Significance level (default = 0.05)
#' 
#' @return A data frame containing:
#'   - theta0: True proportion values
#'   - type1_error: Estimated Type I error rate
#'   - mean_time_ms: Average computation time in milliseconds
#' 
#' @details
#' For each θ₀ in the input vector:
#' 1. Generates B binomial samples from Binom(n, θ₀)
#' 2. Performs exact binomial test (one-sided greater than)
#' 3. Computes proportion of p-values < α (Type I error)
#' 4. Tracks computation times
#' 
#' Includes progress bar and interim results messages.
compute_exact_test_properties <- function(theta0, n, B, alpha = 0.05) {
  # Initialize results data frame
  results <- data.frame(
    theta0 = theta0,
    type1_error = numeric(length(theta0)),
    mean_time_ms = numeric(length(theta0))
  )
  
  # Loop through theta0 values
  for (j in seq_along(theta0)) {
    p_values <- numeric(B)  # Store p-values
    compute_times <- numeric(B)  # Store computation times
    
    # Initialize progress bar
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    
    # Monte Carlo replications
    for (i in 1:B) {
      set.seed(i)  # Reproducibility
      start_time <- Sys.time()
      
      # Generate data under null
      X <- rbinom(1, size = n, prob = theta0[j])
      
      # Perform exact binomial test (one-sided greater than)
      test_result <- binom.test(X, n, p = theta0[j], alternative = "greater")
      p_values[i] <- test_result$p.value
      
      # Record computation time (milliseconds)
      compute_times[i] <- as.numeric(Sys.time() - start_time) * 1000
      
      # Update progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)  # Close progress bar
    
    # Calculate Type I error and average time
    results$type1_error[j] <- mean(p_values < alpha)
    results$mean_time_ms[j] <- mean(compute_times)
    
    # Print interim results
    message(sprintf("theta0 = %.1f | Type I Error = %.3f | Mean Time = %.2f ms",
                    theta0[j], results$type1_error[j], results$mean_time_ms[j]))
  }
  
  return(results)
}

# Simulation parameters
theta0 <- seq(0.1, 0.9, by = 0.1)  # True proportions to evaluate
n <- 100                            # Fixed sample size
B <- 10^4                            # Number of replications
alpha <- 0.05                       # Significance level

# Run simulation
results <- compute_exact_test_properties(theta0, n, B, alpha)

# View results
print(results)

# Plot Type I error rates
plot(results$theta0, results$type1_error, type = "b",
     xlab = "True Proportion (θ₀)", ylab = "Type I Error Rate",
     main = "Exact Binomial Test Type I Error Rates",
     ylim = c(0, max(results$type1_error)*1.1))
abline(h = alpha, col = "red", lty = 2)

# Optional: Save results
# write.csv(results, "exact_test_properties.csv", row.names = FALSE)
