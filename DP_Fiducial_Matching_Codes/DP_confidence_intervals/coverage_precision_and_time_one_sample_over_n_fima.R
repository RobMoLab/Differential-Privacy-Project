##############################################################################
# Fiducial Approach for Differentially Private Proportion Estimation
#
# Description: Implements and evaluates a fiducial-based approach for computing
#              differentially private confidence intervals for binomial proportions
#
# Functions:
#   get_pi - Computes private estimate of proportion
#   fima - Fiducial interval Monte Carlo approximation
#   One_sample_CI_Fiducial_Appl - Computes CI for one-sample proportion
#   precision_coverage_time - Evaluates performance across sample sizes
#
# Dependencies: ggplot2, tidyr
#
# Author: [Your Name]
# Date: [Date]
# Version: 1.0
###############################################################################

library(ggplot2)
library(tidyr)

get_pi <- function(theta, u, w, eps, n) {
  # Input validation
  if (theta <= 0 || theta >= 1) stop("theta must be between 0 and 1")
  if (eps <= 0) stop("eps must be positive")
  if (n <= 0) stop("n must be positive")
  if (any(u < 0 | u > 1)) stop("u must be between 0 and 1")
  if (abs(w) > 0.5) stop("w must be between -0.5 and 0.5")
  
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}

#Fiducial Interval Monte Carlo Approximation

fima <- function(pi0, eps, n, H = 1, delta = 1, seed = 123) {
  # Input validation
  if (eps <= 0) stop("eps must be positive")
  if (n <= 0) stop("n must be positive")
  if (H <= 0) stop("H must be positive")
  
  set.seed(seed)
  Wj <- runif(H, -0.5, 0.5)
  pi_star <- pi0 + delta/(eps*n)*sign(Wj)*log(1 - 2*abs(Wj))
  
  theta_ifb <- rep(NA, H)
  index <- pi_star < 1 & pi_star > 0
  
  tryCatch({
    theta_ifb[index] <- rbeta(
      sum(index), 
      n*pi_star[index] + 0.5, 
      n*(1 - pi_star[index]) + 0.5
    )
  }, error = function(e) {
    stop("Error in beta sampling: ", e$message)
  })
  
  theta_ifb[pi_star >= 1] <- 1 - .Machine$double.eps
  theta_ifb[pi_star <= 0] <- .Machine$double.eps
  
  return(theta_ifb)
}

#Compute CI for One-sample Proportion Using Fiducial Approach

One_sample_CI_Fiducial_Appl <- function(theta0=0.7, H=100, eps=1, n=30, alpha=0.05, seed=123) {
  # Input validation
  if (theta0 <= 0 || theta0 >= 1) stop("theta0 must be between 0 and 1")
  if (H <= 0) stop("H must be positive")
  if (eps <= 0) stop("eps must be positive")
  if (n <= 0) stop("n must be positive")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  
  tryCatch({
    pi0 <- get_pi(theta0, u=runif(n), w=runif(1, -0.5, 0.5), eps=eps, n=n)
    JINI_solution <- fima(pi0=pi0, H=H, eps=eps, n=n, seed=seed)
    c(quantile(JINI_solution, alpha/2, na.rm=TRUE), 
      quantile(JINI_solution, 1-(alpha/2), na.rm=TRUE))
  }, error = function(e) {
    warning("Error in CI computation: ", e$message)
    return(c(NA, NA))
  })
}

#Evaluate Precision and Coverage Across Sample Sizes

precision_coverage_time <- function(H, theta0, nvalues, eps, seed, B, alpha) {
  # Input validation
  if (H <= 0) stop("H must be positive")
  if (theta0 <= 0 || theta0 >= 1) stop("theta0 must be between 0 and 1")
  if (any(nvalues <= 0)) stop("All sample sizes must be positive")
  if (eps <= 0) stop("eps must be positive")
  if (B <= 0) stop("B must be positive")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  
  # Pre-allocate results list for better performance
  results <- vector("list", length(nvalues))
  
  for (i in seq_along(nvalues)) {
    n <- nvalues[i]
    cat("Processing n =", n, "\n")
    
    # Initialize storage
    times <- numeric(B)
    lengths <- numeric(B)
    covers <- logical(B)
    
    for (j in 1:B) {
      tryCatch({
        set.seed(seed + j)
        
        start_time <- Sys.time()
        CI <- One_sample_CI_Fiducial_Appl(
          theta0=theta0, H=H, eps=eps, n=n, 
          alpha=alpha, seed=seed+j
        )
        end_time <- Sys.time()
        
        times[j] <- as.numeric(end_time - start_time) * 1000
        lengths[j] <- CI[2] - CI[1]
        covers[j] <- (CI[1] <= theta0) & (CI[2] >= theta0)
      }, error = function(e) {
        message(sprintf("Error in replication %d for n=%d: %s", j, n, e$message))
        times[j] <- NA
        lengths[j] <- NA
        covers[j] <- NA
      })
    }
    
    # Store results for this sample size
    results[[i]] <- data.frame(
      n = n,
      mean_time = mean(times, na.rm=TRUE),
      median_time = median(times, na.rm=TRUE),
      mean_length = mean(lengths, na.rm=TRUE),
      coverage = mean(covers, na.rm=TRUE),
      n_valid = sum(!is.na(covers)),
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results and add metadata
  results <- do.call(rbind, results)
  attr(results, "simulation_params") <- list(
    theta0 = theta0,
    eps = eps,
    alpha = alpha,
    H = H,
    B = B,
    seed = seed,
    timestamp = Sys.time(),
    R_version = R.version.string
  )
  return(results)
}

# Run the evaluation with timing
sim_time <- system.time({
  Coverage_Precision_and_time_one_sample_over_n_values <- precision_coverage_time(
    H = 10^3, 
    theta0 = 0.2, 
    nvalues = c(8*2^(0:6), 700, 900, 1200), 
    eps = 1, 
    seed = 123, 
    B = 10^4, 
    alpha = 0.05
  )
})

# Print summary information
cat("\nSimulation completed in", round(sim_time["elapsed"], 1), "seconds\n")
print(attr(Coverage_Precision_and_time_one_sample_over_n_values, "simulation_params"))
print(Coverage_Precision_and_time_one_sample_over_n_values)

