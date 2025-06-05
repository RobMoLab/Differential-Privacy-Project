###############################################################################
# Differentially Private Confidence Interval Simulation
#
# Description: Evaluates coverage and properties of differentially private
#              binomial CIs using Tulap noise across varying sample sizes.
#
# Functions:
#   tulap_ci - Computes DP CI using Tulap mechanism
#   sdp_fun - Helper function for binomial SDP mechanism
#   compute_coverage_varying_n - Runs simulation study across sample sizes
#
# Dependencies: binomialDP

###############################################################################

library(binomialDP)


# Compute Differentially Private CI using Tulap Mechanism

tulap_ci <- function(theta0, n, eps, alpha) {
  # Input validation
  if (theta0 <= 0 || theta0 >= 1) stop("theta0 must be between 0 and 1")
  if (n <= 0) stop("n must be positive")
  if (eps <= 0) stop("eps must be positive")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  
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
  # Input validation
  if (seed < 0 || seed > 1) stop("seed must be between 0 and 1")
  if (theta <= 0 || theta >= 1) stop("theta must be between 0 and 1")
  if (n <= 0) stop("n must be positive")
  
  B <- qbinom(seed, size = n, prob = theta)
  s1 <- B + (1/ep) * N[1]
  return(s1)
}

# Run Simulation Study Across Sample Sizes for DP CIs
compute_coverage_varying_n <- function(theta_value, sample_sizes, eps, alpha, n_replicates) {
  # Input validation
  if (theta_value <= 0 || theta_value >= 1) stop("theta_value must be between 0 and 1")
  if (any(sample_sizes <= 0)) stop("All sample sizes must be positive")
  if (eps <= 0) stop("eps must be positive")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  if (n_replicates <= 0) stop("n_replicates must be positive")
  
  # Pre-allocate results list for better performance
  results <- vector("list", length(sample_sizes))
  
  for (j in seq_along(sample_sizes)) {
    n <- sample_sizes[j]
    coverage_count <- 0
    ci_lengths <- numeric(n_replicates)
    run_times <- numeric(n_replicates)
    valid_count <- 0
    
    pb <- utils::txtProgressBar(min = 0, max = n_replicates, style = 3)
    
    for (i in 1:n_replicates) {
      tryCatch({
        set.seed(i)  # For reproducibility
        
        # Time the CI computation
        start_time <- Sys.time()
        ci <- tulap_ci(theta_value, n, eps, alpha)
        end_time <- Sys.time()
        
        # Store results
        run_times[i] <- as.numeric(end_time - start_time) * 1000
        ci_lengths[i] <- ci[2] - ci[1]
        
        if (ci[1] <= theta_value && theta_value <= ci[2]) {
          coverage_count <- coverage_count + 1
        }
        valid_count <- valid_count + 1
      }, error = function(e) {
        ci_lengths[i] <- NA
        run_times[i] <- NA
        warning(sprintf("Error in replication %d for n=%d: %s", i, n, e$message))
      })
      utils::setTxtProgressBar(pb, i)
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
  
  # Combine results and add metadata
  results <- do.call(rbind, results)
  attr(results, "simulation_params") <- list(
    theta_value = theta_value,
    eps = eps,
    alpha = alpha,
    n_replicates = n_replicates,
    timestamp = Sys.time(),
    package_versions = list(
      binomialDP = packageVersion("binomialDP"),
      R = R.version.string
    )
  )
  return(results)
}

# Run the simulation with timing
theta_value <- 0.2
sample_sizes <- c(8 * 2^(0:6), 700, 900, 1200)
eps <- 1
alpha <- 0.05
n_replicates <- 10^4  # Reduced for testing

sim_time <- system.time({
  results <- compute_coverage_varying_n(
    theta_value = theta_value,
    sample_sizes = sample_sizes,
    eps = eps,
    alpha = alpha,
    n_replicates = n_replicates
  )
})

# Print simulation summary
cat("\nSimulation completed in", round(sim_time["elapsed"], 1), "seconds\n")
print(attr(results, "simulation_params"))
print(results)


