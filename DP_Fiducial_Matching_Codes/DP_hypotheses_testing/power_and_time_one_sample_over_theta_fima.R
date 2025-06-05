###############################################################################
# DIFFERENTIALLY PRIVATE HYPOTHESIS TESTING FOR BINOMIAL PROPORTIONS
#
# Description: Implements and evaluates a differentially private hypothesis test
#              for binomial proportions using a fiducial approach
#
# Functions:
#   get_pi - Computes private estimate of proportion
#   fima - Fiducial interval Monte Carlo approximation
#   compute_type1_error - Evaluates type I error rate across parameter space
###############################################################################

#' Compute Differentially Private Estimate of Proportion

get_pi <- function(theta, u, w, eps, n) {
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}

#' FIMA algorithm
#'
#' Generates samples from the fiducial distribution for differentially private
#' binomial proportion estimation.
#'
#' @param pi0 Initial private estimate
#' @param eps Privacy parameter (epsilon > 0)
#' @param n Sample size (positive integer)
#' @param H Number of Monte Carlo samples (positive integer, default=1)
#' @param seed Random seed (default=123)
#' @return Vector of H fiducial samples

fima <- function(pi0, eps, n, H = 1, delta = 1, seed = 123) {
  set.seed(seed)
  
  Wj <- runif(H, -0.5, 0.5)
  pi_star <- pi0 + delta/(eps*n)*sign(Wj)*log(1 - 2*abs(Wj))
  
  theta_ifb <- rep(NA, H)
  
  index <- pi_star < 1 & pi_star > 0
  theta_ifb[index] <- rbeta(length(index), 
                            n*pi_star[index] + 0.5, 
                            n*(1 - pi_star[index]) + 0.5)
  theta_ifb[pi_star >= 1] <- 1 - .Machine$double.eps
  theta_ifb[pi_star <= 0] <- .Machine$double.eps
  
  return(theta_ifb)
}

###############################################################################
# SIMULATION PARAMETERS DOCUMENTATION
#
# Parameters for Type I Error Evaluation:
#   H = 10^3 - Number of fiducial samples per test
#
#   B = 10^4 - Number of test replications
#
#   n = 100 - Sample size
#
#   theta = seq(0.2,0.9,0.1) - True proportion values to evaluate
#
#   eps = 1 - Privacy budget
#
#   alpha = 0.05 - Significance level
#
#   theta0 = 0.2 - Null hypothesis value (H0: θ ≤ θ0)

###############################################################################

#' Evaluate power of DP Hypothesis Test
#'
#' Performs simulation study to evaluate the power of a
#' differentially private hypothesis test for binomial proportions using Fima.
#'
#' @param theta0 Null hypothesis value (H0: θ ≤ θ0)
#' @param theta Vector of true proportion values to evaluate
#' @param n Sample size (positive integer)
#' @param eps Privacy parameter (epsilon > 0)
#' @param B Number of test replications (positive integer)
#' @param H Number of fiducial samples per test (positive integer)
#' @return Data frame with columns:
#'   - theta: True proportion value
#'   - power_fima: power of test using Fima (power when θ = θ0)
#'   - mean_time_fima: Mean computation time per test (ms)
#'   - median_time_fima: Median computation time per test (ms)
#'

compute_type1_error <- function(theta0, theta, n, eps, B, H) {
  results <- data.frame(
    theta = numeric(),
    power_fima = numeric(),
    mean_time_fima = numeric(),
    median_time_fima = numeric(),
    stringsAsFactors = FALSE
  )
  
  res_jini <- matrix(NA, B, length(theta))
  compute_times <- matrix(NA, B, length(theta))
  
  for (j in 1:length(theta)) {
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    
    for (i in 1:B) {
      set.seed(i)
      
      # Time the computation
      start_time <- Sys.time()
      
      pi0 <- get_pi(theta = theta[j], u = runif(n), 
                    w = runif(1, -0.5, 0.5), eps = eps, n = n)
      
      # JINI P-Value
      emp_density <- fima(pi0 = pi0, H = H, eps = eps, 
                          n = n, seed = i + 2*H)
      res_jini[i,j] <- sum(emp_density <= theta0)/(H + 1)
      
      end_time <- Sys.time()
      compute_times[i,j] <- as.numeric(end_time - start_time)*1000 #covert time to ms
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Calculate the power and timing statistics
    results <- rbind(results, data.frame(
      theta = theta[j],
      power_fima = mean(res_jini[,j] < alpha),
      mean_time_fima = mean(compute_times[,j]),
      median_time_fima = median(compute_times[,j])
    ))
    
    print(paste("Completed theta =", theta[j]))
  }
  
  return(results)
}

# Simulation parameters
H <- 10^3  # Number of fiducial samples per test
B <- 10^4  # Number of test replications  
n <- 100   # Sample size
theta <- seq(0.2, 0.9, 0.1)  # True proportion values to evaluate
eps <- 1    # Privacy budget
alpha <- 0.05  # Significance level
theta0 <- 0.2  # Null hypothesis value

# Run simulation
results <- compute_type1_error(theta0, theta, n, eps, B, H)

# View results
print(results)
