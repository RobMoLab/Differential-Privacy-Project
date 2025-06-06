#' Differentially Private Proportion Estimator
#'
#' Computes a private estimate of a proportion using the Laplace mechanism.
#' Implements ε-differential privacy.
#'
#' @param theta True proportion
#' @param u Vector of uniform random variables (length n)
#' @param w Uniform random variable in [-0.5, 0.5] for noise generation
#' @param eps Privacy parameter (ε)
#' @param n Sample size
#'
#' @return Differentially private estimate of the proportion
#' @examples
#' get_pi(theta = 0.5, u = runif(100), w = runif(1, -0.5, 0.5), eps = 1, n = 100)
get_pi <- function(theta, u, w, eps, n) {
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}



#'FIMA
#'
#' Generates fiducial samples for differentially private inference using
#' a fiducial approach with beta distributions.
#'
#' @param pi0 Private proportion estimate
#' @param eps Privacy parameter (ε)
#' @param n Sample size
#' @param H Number of fiducial samples (default = 1)
#' @param delta Sensitivity parameter (default = 1)
#' @param seed Random seed (default = 123)
#'
#' @return Vector of H fiducial samples
#' @examples
#' fima(pi0 = 0.5, eps = 1, n = 100, H = 1000)
fima <- function(pi0, eps, n, H = 1, delta = 1, seed = 123) {
  
  set.seed(seed)
  
  # Generate noise variables
  Wj <- runif(H, -0.5, 0.5)
  
  # Create perturbed proportions
  pi_star <- pi0 + delta/(eps*n)*sign(Wj)*log(1 - 2*abs(Wj))
  
  # Initialize fiducial samples
  theta_ifb <- rep(NA, H)
  
  # Generate beta samples for valid proportions
  index <- pi_star < 1 & pi_star > 0
  theta_ifb[index] <- rbeta(
    sum(index), 
    n*pi_star[index] + 0.5,  # Add 0.5 for Haldane-Anscombe correction
    n*(1 - pi_star[index]) + 0.5
  )
  
  # Handle boundary cases
  theta_ifb[pi_star >= 1] <- 1 - .Machine$double.eps
  theta_ifb[pi_star <= 0] <- .Machine$double.eps
  
  return(theta_ifb)
}



#' Compute Type I Error for DP Fiducial Test
#'
#' Evaluates Type I error rates of the differentially private fiducial test
#' across different null proportions θ₀.
#'
#' @param theta0 Vector of null proportions to evaluate
#' @param n Sample size
#' @param eps Privacy parameter (ε)
#' @param B Number of Monte Carlo replications
#' @param H Number of fiducial samples
#'
#' @return Data frame with:
#'   - theta: Null proportion values
#'   - type1_error: Estimated Type I error rate
#'   - mean_time: Average computation time (ms)
#'   - median_time: Median computation time (ms)
#' @examples
#' compute_type1_error(theta0 = seq(0.1, 0.9, by = 0.1), n = 100, eps = 1, B = 100, H = 1000)
compute_type1_error <- function(theta0, n, eps, B, H) {
  # Initialize results structure
  results <- data.frame(
    theta = numeric(),
    type1_error = numeric(),
    mean_time = numeric(),
    median_time = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Initialize storage matrices
  res_jini <- matrix(NA, B, length(theta0))
  compute_times <- matrix(NA, B, length(theta0))
  
  # Main simulation loop
  for (j in 1:length(theta0)) {
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    
    for (i in 1:B) {
      set.seed(i)
      start_time <- Sys.time()
      
      # Generate private estimate
      pi0 <- get_pi(
        theta = theta0[j], 
        u = runif(n), 
        w = runif(1, -0.5, 0.5),
        eps = eps,
        n = n
      )
      
      # Compute fiducial p-value
      emp_density <- fima(
        pi0 = pi0,
        H = H,
        eps = eps,
        n = n,
        seed = i + 2*H
      )
      res_jini[i,j] <- sum(emp_density > theta0[j])/(H + 1)
      
      # Record computation time (ms)
      end_time <- Sys.time()
      compute_times[i,j] <- as.numeric(end_time - start_time)*1000
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Store results for this theta0
    results <- rbind(results, data.frame(
      theta = theta0[j],
      type1_error = mean(res_jini[,j] < 0.05),
      mean_time = mean(compute_times[,j]),
      median_time = median(compute_times[,j])
    ))
    
    print(paste("Completed theta =", theta0[j]))
  }
  
  return(results)
}



# Simulation parameters
H <- 10^3       # Number of fiducial samples
B <- 10^4       # Number of Monte Carlo replications
n <- 100        # Sample size
theta0 <- seq(0.1, 0.9, by = 0.1)  # Null proportions to evaluate
eps <- 1        # Privacy parameter (ε = 1)

# Run simulation
results <- compute_type1_error(theta0, n, eps, B, H)

# View results
print(results)

# Plot Type I error rates
plot(results$theta, results$type1_error, type = "b",
     xlab = "Null Proportion (θ₀)", ylab = "Type I Error Rate",
     main = "DP Fiducial Test Type I Error Rates",
     ylim = c(0, max(results$type1_error)*1.1))
abline(h = 0.05, col = "red", lty = 2)

# Save results
#write.csv(results, "dp_fiducial_type1_error.csv", row.names = FALSE)
