

# Compute Non-Private Wald Confidence Interval for Binomial Proportion

Non_private_CI <- function(p, n, alpha) {
  # Simulate Bernoulli trials
  observed_sample <- rbinom(n, 1, p)
  success_proportion <- sum(observed_sample) / n
  
  # Compute confidence interval
  se <- sqrt(success_proportion * (1 - success_proportion) / n)
  margin_of_error <- qnorm(1 - alpha / 2) * se
  ci_lower <- success_proportion - margin_of_error
  ci_upper <- success_proportion + margin_of_error
  
  return(c(ci_lower, ci_upper))
}

#' Evaluate Coverage Properties of Wald Confidence Intervals
#'
#' Performs simulation study to evaluate coverage probability, CI lengths,
#' and computation time of Wald CIs across a range of true proportions.
#'
#' @param thetavalues Vector of true proportion values to evaluate
#' @param n Sample size (positive integer)
#' @param seed Random seed for reproducibility
#' @param B Number of replications per theta value
#' @param alpha Significance level
#' @return Data frame with columns:
#'   - theta: True proportion value
#'   - coverage_norm_apprx: Empirical coverage probability
#'   - avg_ci_length_norm_apprx: Average CI width
#'   - mean_time_norm_apprx: Mean computation time (ms)
#'
#' @details
#' For each theta in thetavalues:
#' 1. Repeats B times:
#'    a. Simulates binomial data and computes Wald CI
#'    b. Records whether CI covers true theta
#'    c. Records CI width and computation time
#' 2. Calculates summary statistics across replications
#' 3. Returns results in tidy data frame
#'
#' Includes progress bar to track simulation progress.
CI_coverage_comparison <- function(thetavalues, n, seed, B, alpha) {
  results <- data.frame(
    theta = numeric(),
    coverage = numeric(),
    avg_ci_length = numeric(),
    mean_time = numeric(),
    median_time = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (theta in thetavalues) {
    coverage_count <- 0
    ci_lengths <- numeric(B)
    run_times <- numeric(B)
    
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    
    for (h in 1:B) {
      set.seed(seed + h)
      
      # Time the CI computation
      start_time <- Sys.time()
      CI_NP <- Non_private_CI(p = theta, n = n, alpha = alpha)
      end_time <- Sys.time()
      
      # Store results
      run_times[h] <- as.numeric(end_time - start_time)*1000
      ci_lengths[h] <- CI_NP[2] - CI_NP[1]
      
      # Check coverage
      if (CI_NP[1] <= theta && theta <= CI_NP[2]) {
        coverage_count <- coverage_count + 1
      }
      setTxtProgressBar(pb, h)
    }
    close(pb)
    
    # Calculate summary statistics
    coverage_prob <- coverage_count / B
    time_stats <- c(mean = mean(run_times), 
                    median = median(run_times),
                    sd = sd(run_times))
    length_stats <- c(mean = mean(ci_lengths), 
                      median = median(ci_lengths),
                      sd = sd(ci_lengths))
    
    # Store results for this theta
    results <- rbind(results, data.frame(
      theta = theta,
      coverage_norm_apprx = coverage_prob,
      avg_ci_length_norm_apprx= length_stats["mean"],
      mean_time_norm_apprx = time_stats["mean"]
    ))
  }
  
  return(results)
}

# Simulation parameters
theta_values <- seq(from = 0.1, to = 0.985, by = 0.01)
n <- 30
alpha <- 0.05
B <- 10^4  # Number of replications
seed <- 123

# Run simulation
results <- CI_coverage_comparison(
  thetavalues = theta_values,
  n = n,
  seed = seed,
  B = B,
  alpha = alpha
)

# View results
results
