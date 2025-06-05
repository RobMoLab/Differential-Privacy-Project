###############################################################################
# Exact Binomial Test Coverage Simulation
#
# Description: This script simulates the coverage probability and confidence 
# interval properties of exact binomial tests across different true probability 
# values (theta). It evaluates how often the exact binomial confidence intervals
# contain the true probability parameter and examines the interval widths and
# computation times.
#
# Inputs:
#   - theta_values: Vector of true probability values to evaluate (0-1)
#   - n: Number of Bernoulli trials (positive integer)
#   - alpha: Significance level (0-1), typically 0.05
#   - n_replicates: Number of simulation replicates per theta value (positive integer)
#
# Outputs:
#   Data frame containing for each theta:
#   - theta: Input probability value
#   - coverage_exac_binom: Empirical coverage probability
#   - avg_ci_length_exac_binom: Average confidence interval width
#   - mean_time_exac_binom: Average computation time in milliseconds
#
######################################################################################################################


###################### Exact Binomial test Coverage and precision, while varying theta ##########
#function to compute coverage with timing
simulate_binom_coverage <- function(theta_values, n, alpha, n_replicates) {
  results <- data.frame(
    theta = numeric(),
    coverage = numeric(),
    avg_ci_length = numeric(),
    mean_time = numeric(),
    median_time = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (theta in theta_values) {
    coverage_count <- 0
    ci_lengths <- numeric(n_replicates)
    run_times <- numeric(n_replicates)
    
    pb <- txtProgressBar(min = 0, max = n_replicates, style = 3)
    
    for (i in 1:n_replicates) {
      set.seed(i)  # For reproducibility
      
      # Time the CI computation
      start_time <- Sys.time()
      x <- rbinom(1, size = n, prob = theta)
      test <- binom.test(x, n, conf.level = 1 - alpha)
      ci <- test$conf.int
      end_time <- Sys.time()
      
      # Store results
      run_times[i] <- as.numeric(end_time - start_time)*1000  # convert to milliseconds
      ci_lengths[i] <- ci[2] - ci[1]
      
      # Check coverage
      if (ci[1] <= theta && theta <= ci[2]) {
        coverage_count <- coverage_count + 1
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Calculate summary statistics
    coverage_prob <- coverage_count / n_replicates
    time_stats <- c(mean = mean(run_times), 
                    median = median(run_times),
                    sd = sd(run_times))
    length_stats <- c(mean = mean(ci_lengths), 
                      median = median(ci_lengths),
                      sd = sd(ci_lengths))
    
    # Store results for this theta
    results <- rbind(results, data.frame(
      theta = theta,
      coverage_exac_binom = coverage_prob,
      avg_ci_length_exac_binom = length_stats["mean"],
      mean_time_exac_binom = time_stats["mean"]
      #median_time = time_stats["median"]
    ))
  }
  
  return(results)
}

# Run simulation
theta_values <- seq(from = 0.1, to = 0.985, by = 0.01)
n <- 30
alpha <- 0.05
n_replicates <- 10^4  

results <- simulate_binom_coverage(theta_values, n, alpha, n_replicates)
