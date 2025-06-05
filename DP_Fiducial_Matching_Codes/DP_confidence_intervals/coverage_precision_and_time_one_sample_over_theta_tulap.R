###############################################################################
# DIFFERENTIALLY PRIVATE CONFIDENCE INTERVALS USING TULAP NOISE
#
# Description: Implements and evaluates differentially private confidence intervals
#              for binomial proportions using Tulap noise mechanism
#
# Dependencies: binomialDP
###############################################################################

library(binomialDP)


###############################################################################
# SIMULATION PARAMETERS DOCUMENTATION
#
# Parameters:
#   theta_values - Sequence of true proportion values (0.1 to 0.985 by 0.01)
#     
#
#   n = 30 - Sample size for binomial data
#     
#
#   eps = 1 - Privacy budget (epsilon)
#     
#
#   alpha = 0.05 - Significance level
#     
#
#   n_replicates = 10^4 - Number of simulation replications

###############################################################################

# Compute Differentially Private CI Using Tulap Mechanism

tulap_ci <- function(theta0, n, eps, alpha) {
  de <- 0.0                    # Delta parameter for Tulap
  b <- exp(-eps)               # Tulap parameter
  q <- 2 * de * b / (1 - b + 2 * de * b)  # Tulap parameter
  seed <- runif(1, 0, 1)
  N <- rtulap(n = 1, m = 0, b = b, q = q)
  Z <- sdp_fun(seed, N, eps, theta0, n)
  ci <- CITwoSide(alpha = alpha, Z, size = n, b = b, q = q)
  return(c(ci[1], ci[2]))
}

# Helper Function for Binomial SDP Mechanism

sdp_fun <- function(seed, N, ep, theta, n) {
  B <- qbinom(seed, size = n, prob = theta)
  s1 <- B + (1/ep) * N[1]
  return(s1)
}

#' Evaluate Coverage Properties of DP Confidence Intervals
#'
#' Performs simulation study to evaluate coverage probability, CI lengths,
#' and computation time of DP CIs across a range of true proportions.
#'
#' @param theta_values Vector of true proportion values to evaluate
#' @param n Sample size (positive integer)
#' @param eps Privacy parameter (epsilon > 0)
#' @param alpha Significance level (0 < alpha < 1)
#' @param n_replicates Number of simulation replications (positive integer)
#' @return Data frame with columns:
#'   - theta: True proportion value
#'   - coverage_tulap: Empirical coverage probability
#'   - mean_time_tulap: Mean computation time (ms)
#'   - mean_length_tulap: Average CI width
#'
#' @details
#' For each theta in theta_values:
#' 1. Repeats n_replicates times:
#'    a. Computes DP CI using tulap_ci()
#'    b. Records whether CI covers true theta
#'    c. Records CI width and computation time
#' 2. Calculates summary statistics across replications
#' 3. Returns results in tidy data frame
compute_coverage <- function(theta_values, n, eps, alpha, n_replicates = 10000) {
  results <- data.frame()
  timing_results <- list()
  
  for (theta in theta_values) {
    coverage_count <- 0
    ci_lengths <- numeric(n_replicates)
    run_times <- numeric(n_replicates)
    
    for (i in 1:n_replicates) {
      set.seed(i)  # For reproducibility
      
      # Time the CI computation
      start_time <- Sys.time()
      ci <- tulap_ci(theta, n, eps, alpha)
      end_time <- Sys.time()
      
      # Store results
      run_times[i] <- as.numeric(end_time - start_time)*1000
      ci_lengths[i] <- ci[2] - ci[1]
      
      # Check coverage
      if (ci[1] <= theta && theta <= ci[2]) {
        coverage_count <- coverage_count + 1
      }
    }
    
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
      coverage_tulap = coverage_prob,
      mean_time_tulap = time_stats["mean"],
      mean_length_tulap = length_stats["mean"]
    ))
    
    timing_results[[as.character(theta)]] <- run_times
  }
  
  return(results)
}

# Simulation parameters
theta_values <- seq(from = 0.1, to = 0.985, by = 0.01)
n <- 30
eps <- 1
alpha <- 0.05
n_replicates <- 10^4 

# Run simulation
results <- compute_coverage(theta_values, n, eps, alpha, n_replicates)

# View results
results
