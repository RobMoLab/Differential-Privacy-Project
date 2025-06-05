###############################################################################
# Non-Private Confidence Interval Simulation
#
# Description: Evaluates coverage and properties of standard binomial CIs
#              across varying sample sizes.
#
# Functions:
#   Non_private_CI - Computes Wald CI for binomial proportion
#   compute_metrics - Runs simulation study across sample sizes

###############################################################################

# Non-private CI function with timing
Non_private_CI <- function(p, n, alpha) {
  # Input validation
  if (p <= 0 || p >= 1) stop("p must be between 0 and 1")
  if (n <= 0) stop("n must be positive")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  
  # Simulate Bernoulli trials
  observed_sample <- rbinom(n, 1, p)
  success_proportion <- sum(observed_sample) / n
  
  # Compute confidence interval
  se <- sqrt(success_proportion * (1 - success_proportion) / n)
  margin_of_error <- qnorm(1 - alpha / 2) * se
  ci_lower <- success_proportion - margin_of_error
  ci_upper <- success_proportion + margin_of_error
  
  return(list(ci = c(ci_lower, ci_upper), 
              success_proportion = success_proportion))
}

# Function to compute all metrics
compute_metrics <- function(theta0, nvalues, B, alpha) {
  # Pre-allocate results list for better performance
  results <- vector("list", length(nvalues))
  
  for (j in seq_along(nvalues)) {
    n <- nvalues[j]
    coverage_count <- 0
    ci_lengths <- numeric(B)
    run_times <- numeric(B)
    
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    
    for (i in 1:B) {
      tryCatch({
        set.seed(i)
        
        # Time the CI computation
        start_time <- Sys.time()
        res <- Non_private_CI(p = theta0, n = n, alpha = alpha)
        end_time <- Sys.time()
        
        # Store results
        run_times[i] <- as.numeric(end_time - start_time)*1000 # convert to milliseconds
        ci_lengths[i] <- res$ci[2] - res$ci[1]
        
        # Check coverage
        if (res$ci[1] <= theta0 && theta0 <= res$ci[2]) {
          coverage_count <- coverage_count + 1
        }
      }, error = function(e) {
        ci_lengths[i] <- NA
        run_times[i] <- NA
      })
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Calculate summary statistics
    results[[j]] <- data.frame(
      n = n,
      coverage = coverage_count / B,
      avg_ci_length = mean(ci_lengths, na.rm = TRUE),
      mean_time = mean(run_times, na.rm = TRUE),
      median_time = median(run_times, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all results and add metadata
  results <- do.call(rbind, results)
  attr(results, "simulation_params") <- list(
    theta0 = theta0,
    alpha = alpha,
    replications = B,
    timestamp = Sys.time()
  )
  return(results)
}

# Run simulation
theta0 <- 0.2   # Fixed parameter of interest
nvalues <- c(8*2^(0:6), 700, 900, 1200)   # sample sizes
B <- 10^4  # Number of replications
alpha <- 0.05

results <- compute_metrics(
  theta0 = theta0,
  nvalues = nvalues,
  B = B,
  alpha = alpha
)

# Print simulation parameters
print(attr(results, "simulation_params"))
print(results)
