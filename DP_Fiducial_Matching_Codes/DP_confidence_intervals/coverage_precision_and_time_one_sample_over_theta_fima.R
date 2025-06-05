###############################################################################
# Differentially Private Fiducial Confidence Intervals for Binomial Proportions
#
# Description: Implements and evaluates a fiducial-based approach for computing
#              differentially private confidence intervals for binomial proportions
#              across a range of true proportion values.
#
# Functions:
#   get_pi - Computes private estimate of proportion
#   fima - Fiducial interval Monte Carlo approximation
#   One_sample_CI_Fiducial - Computes CI for one-sample proportion
#   CI_coverage_comparison - Evaluates coverage and precision across theta values
#
# Dependencies:  progress (optional)

###############################################################################


# Compute Private Estimate of Proportion

get_pi <- function(theta, eps, n) {
  # Input validation
  if (theta <= 0 || theta >= 1) stop("theta must be between 0 and 1")
  if (eps <= 0) stop("eps must be positive")
  if (n <= 0) stop("n must be positive")
  
  w <- runif(1, -0.5, 0.5)
  u <- runif(n)
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}

# FIMA algorithm 

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
  
  # Handle boundary cases with machine precision
  theta_ifb[pi_star >= 1] <- 1 - .Machine$double.eps
  theta_ifb[pi_star <= 0] <- .Machine$double.eps
  
  return(theta_ifb)
}

# Compute CI for One-sample Proportion Using Fiducial Approach

One_sample_CI_Fiducial <- function(pi0, H, eps, n, alpha, seed) {
  # Input validation
  if (H <= 0) stop("H must be positive")
  if (eps <= 0) stop("eps must be positive")
  if (n <= 0) stop("n must be positive")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  
  tryCatch({
    JINI_solution <- fima(pi0 = pi0, H = H, eps = eps, n = n, seed = seed)
    c(
      quantile(JINI_solution, alpha/2, na.rm = TRUE), 
      quantile(JINI_solution, 1 - (alpha/2), na.rm = TRUE)
    )
  }, error = function(e) {
    warning("Error in CI computation: ", e$message)
    return(c(NA, NA))
  })
}

#' Evaluate Coverage and Precision of Confidence Intervals Across Theta Values
#' 
#' Simulation study to evaluate the performance of FIMA differentially private confidence
#' intervals across a range of true proportion values.
#'
#' @param H Number of fiducial samples
#' @param thetavalues Vector of true proportion values to evaluate
#' @param eps Privacy budget
#' @param n Sample size
#' @param seed Random seed
#' @param B Number of simulation replicates per theta value
#' @param alpha Significance level
#'
#' @return A data frame with simulation results including:
#'   - theta: True proportion value
#'   - coverage_fima: Empirical coverage rate
#'   - mean_length_fima: Average CI width
#'   - median_length_fima: Median CI width  
#'   - mean_time_fima: Average computation time (ms)
#'   - median_time_fima: Median computation time (ms)
#'   - n_valid: Number of successful simulations
#' @export
#'
#' @examples
#' results <- CI_coverage_comparison(
#'   H = 1000,
#'   thetavalues = seq(0.1, 0.9, 0.1),
#'   eps = 1,
#'   n = 30,
#'   seed = 123,
#'   B = 100,
#'   alpha = 0.05
#' )

CI_coverage_comparison <- function(H, thetavalues, eps, n, seed, B, alpha) {
  # Input validation
  if (H <= 0) stop("H must be positive")
  if (any(thetavalues <= 0 | thetavalues >= 1)) stop("All theta values must be between 0 and 1")
  if (eps <= 0) stop("eps must be positive")
  if (n <= 0) stop("n must be positive")
  if (B <= 0) stop("B must be positive")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  
  # Initialize progress bar
  use_progress_pkg <- requireNamespace("progress", quietly = TRUE)
  if (use_progress_pkg) {
    pb <- progress::progress_bar$new(
      format = "  Simulating [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
      total = length(thetavalues) * B,
      clear = FALSE,
      width = 60
    )
  } else {
    pb <- txtProgressBar(min = 0, max = length(thetavalues) * B, style = 3)
    message("Using basic progress bar (install 'progress' package for enhanced bar)")
  }
  
  # Pre-allocate storage
  results <- vector("list", length(thetavalues))
  CI_Lengths <- matrix(NA, nrow = length(thetavalues), ncol = B)
  time_matrix <- matrix(NA, nrow = length(thetavalues), ncol = B)
  
  # Main simulation loop
  for (j in seq_along(thetavalues)) {
    theta <- thetavalues[j]
    coverage_count <- 0
    
    for (h in 1:B) {
      tryCatch({
        set.seed(seed + h)
        start_time <- Sys.time()
        
        pi0 <- get_pi(theta = theta, eps = eps, n = n)
        CI <- One_sample_CI_Fiducial(
          pi0 = pi0, H = H, alpha = alpha, 
          n = n, eps = eps, seed = seed + h
        )
        
        end_time <- Sys.time()
        
        # Store results
        CI_Lengths[j, h] <- CI[2] - CI[1]
        time_matrix[j, h] <- as.numeric(end_time - start_time) * 1000
        coverage_count <- coverage_count + 
          as.numeric(CI[1] <= theta && CI[2] >= theta)
        
        # Update progress
        if (use_progress_pkg) pb$tick() else setTxtProgressBar(pb, (j-1)*B + h)
      }, error = function(e) {
        message(sprintf("Error for theta=%.2f, rep=%d: %s", theta, h, e$message))
      })
    }
    
    # Store summary for this theta
    results[[j]] <- data.frame(
      theta = theta,
      coverage_fima = coverage_count / B,
      mean_length_fima = mean(CI_Lengths[j, ], na.rm = TRUE),
      median_length_fima = median(CI_Lengths[j, ], na.rm = TRUE),
      mean_time_fima = mean(time_matrix[j, ], na.rm = TRUE),
      median_time_fima = median(time_matrix[j, ], na.rm = TRUE),
      n_valid = sum(!is.na(CI_Lengths[j, ])),
      stringsAsFactors = FALSE
    )
  }
  
  # Clean up progress bar
  if (use_progress_pkg) {
    if (!pb$finished) pb$terminate()
  } else {
    close(pb)
  }
  
  # Combine results and add metadata
  results <- do.call(rbind, results)
  attr(results, "simulation_params") <- list(
    H = H,
    eps = eps,
    n = n,
    alpha = alpha,
    B = B,
    seed = seed,
    timestamp = Sys.time(),
    R_version = R.version.string
  )
  
  return(results)
}

# Simulation parameters
thetavalues <- seq(from = 0.1, to = 0.98, by = 0.01)
H <- 10^3
n <- 30
eps <- 1
alpha <- 0.05
B <- 10^4 
seed <- 3

# Run simulation with timing
sim_time <- system.time({
  Coverage_Precision_and_time_one_sample_over_theta_fima <- CI_coverage_comparison(
    H = H, 
    thetavalues = thetavalues, 
    eps = eps, 
    n = n, 
    seed = seed, 
    B = B, 
    alpha = alpha
  )
})

# Print summary information
cat("\nSimulation completed in", round(sim_time["elapsed"], 1), "seconds\n")
#print(attr(Coverage_Precision_and_time_one_sample_over_theta_fima, "simulation_params"))
print(Coverage_Precision_and_time_one_sample_over_theta_fima)

