###############################################################################
# Exact Binomial Test Coverage Simulation
#
# Description: Simulates coverage probability and CI properties for exact binomial
#              tests across varying sample sizes.
#
# Inputs:
#   theta_value  - True probability of success (0 < theta_value < 1)
#   sample_sizes - Vector of sample sizes to evaluate
#   alpha        - Significance level (0 < alpha < 1)
#   n_replicates - Number of simulation replicates per sample size
#
# Output:
#   data.frame with columns:
#     n            - Sample size
#     theta        - True probability
#     coverage     - Empirical coverage probability
#     avg_ci_length- Average confidence interval width
#     mean_time    - Mean computation time (ms)
#     median_time  - Median computation time (ms)
#     n_valid      - Number of valid replicates
#
# Notes:
#   - Uses binom.test() for exact binomial CIs
#   - Includes timing information for performance evaluation
#   - Progress bar shows simulation progress
#
# Example:
#   results <- simulate_binom_coverage_varying_n(
#     theta_value = 0.2,
#     sample_sizes = c(10, 20, 50),
#     alpha = 0.05,
#     n_replicates = 1000
#   )
#
###############################################################################

simulate_binom_coverage_varying_n <- function(theta_value, sample_sizes, alpha, n_replicates) {
  ## Input validation
  if (theta_value <= 0 || theta_value >= 1) stop("theta_value must be between 0 and 1")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  if (any(sample_sizes <= 0)) stop("All sample sizes must be positive")
  
  results <- vector("list", length(sample_sizes))
  
  for (j in seq_along(sample_sizes)) {
    n <- sample_sizes[j]
    coverage_count <- 0
    ci_lengths <- numeric(n_replicates)
    run_times <- numeric(n_replicates)
    valid_count <- 0
    
    pb <- txtProgressBar(min = 0, max = n_replicates, style = 3)
    
    for (i in 1:n_replicates) {
      tryCatch({
        set.seed(i)  # For reproducibility
        
        # Time only the CI computation, not the data generation
        x <- rbinom(1, size = n, prob = theta_value)
        start_time <- Sys.time()
        test <- binom.test(x, n, conf.level = 1 - alpha)
        end_time <- Sys.time()
        
        ci <- test$conf.int
        run_times[i] <- as.numeric(end_time - start_time) * 1000 # convert to milliseconds
        ci_lengths[i] <- ci[2] - ci[1]
        
        if (ci[1] <= theta_value && theta_value <= ci[2]) {
          coverage_count <- coverage_count + 1
        }
        valid_count <- valid_count + 1
      }, error = function(e) {
        ci_lengths[i] <- NA
        run_times[i] <- NA
      })
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Remove NA values from failed runs
    valid_lengths <- ci_lengths[!is.na(ci_lengths)]
    valid_times <- run_times[!is.na(run_times)]
    
    if (valid_count > 0) {
  
      results[[j]] <- data.frame(
        n = n,
        theta = theta_value,
        coverage = coverage_count / valid_count,
        avg_ci_length = mean(valid_lengths),
        mean_time = mean(valid_times),
        median_time = median(valid_times),
        n_valid = valid_count,
        stringsAsFactors = FALSE
      )
    } else {
      warning(sprintf("No valid replications for n = %d", n))
    }
  }
  
  # Combine all results at once
  results <- do.call(rbind, results)
  attr(results, "description") <- "Simulation results for exact binomial test coverage"
  return(results)
}

# Run simulation with parameter validation
theta_value <- 0.2
sample_sizes <- c(8 * 2^(0:6), 700, 900, 1200)
alpha <- 0.05
n_replicates <- 10^4  

## Wrap simulation in system.time for overall timing
sim_time <- system.time({
  results <- simulate_binom_coverage_varying_n(
    theta_value = theta_value,
    sample_sizes = sample_sizes,
    alpha = alpha,
    n_replicates = n_replicates
  )
})

## Print summary of results
cat("\nSimulation completed in", round(sim_time["elapsed"], 1), "seconds\n")
print(results)

