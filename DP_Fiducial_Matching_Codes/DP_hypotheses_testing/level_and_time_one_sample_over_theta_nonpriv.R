# Level Simulation for the test; H0: theta <= theta_0 vs Ha: theta > theta_0

#' Compute Type I Error for Binomial Test with Continuity Correction
#' 
#' Evaluates Type I error rates of a standard binomial test with uniform 
#' continuity correction across different null proportions.
#' 
#' @param theta0 Vector of null proportions to evaluate
#' @param n Sample size (fixed across simulations)
#' @param B Number of Monte Carlo replications per theta0 value
#' 
#' @return A data frame containing:
#'   - theta: Null proportion values
#'   - type1_error: Estimated Type I error rate
#'   - mean_time: Average computation time in seconds
#'   - median_time: Median computation time in seconds
#' 
#' @details
#' For each θ₀ in the input vector:
#' 1. Generates B binomial samples from Binom(n, θ₀)
#' 2. Adds uniform continuity correction (-0.5 to 0.5)
#' 3. Computes exact p-values using uniform CDF
#' 4. Estimates Type I error rate as proportion of p-values < 0.05
#' 5. Tracks computation times
#' 
#' Includes progress bar and completion messages.
compute_type1_error <- function(theta0, n, B) {
  # Initialize results structure
  results <- data.frame(
    theta = numeric(),
    type1_error = numeric(),
    mean_time = numeric(),
    median_time = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Initialize storage matrices
  res_nonpriv <- matrix(NA, B, length(theta0))
  compute_times <- matrix(NA, B, length(theta0))
  
  # Main simulation loop
  for (j in 1:length(theta0)) {
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    
    for (i in 1:B) {
      set.seed(i)  # Ensure reproducibility
      start_time <- Sys.time()
      
      # Generate data under null
      X <- rbinom(n = 1, size = n, prob = theta0[j])
      values <- seq(0, n)
      pdf <- dbinom(values, size = n, prob = theta0[j])
      
      # Add uniform continuity correction
      Z <- X + runif(n = 1, min = -1/2, max = 1/2)
      
      # Compute p-value using uniform CDF
      cdf <- punif(values - Z, min = -1/2, max = 1/2)
      res_nonpriv[i,j] <- t(cdf) %*% pdf
      
      # Record computation time (seconds)
      compute_times[i,j] <- as.numeric(Sys.time() - start_time)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Store results for this theta0
    results <- rbind(results, data.frame(
      theta = theta0[j],
      type1_error = mean(res_nonpriv[,j] < 0.05),
      mean_time = mean(compute_times[,j]),
      median_time = median(compute_times[,j])
    ))
    
    print(paste("Completed theta =", theta0[j]))
  }
  
  return(results)
}


# Simulation parameters
theta0 <- seq(0.1, 0.9, by = 0.1)  # Null proportions to evaluate
n <- 100                            # Fixed sample size
B <- 10^4                           # Number of replications
alpha <- 0.05                       # Significance level

# Run simulation
results <- compute_type1_error(theta0, n, B)

# View results
print(results)


# Optional: Save results
# write.csv(results, "binomial_type1_error_results.csv", row.names = FALSE)
