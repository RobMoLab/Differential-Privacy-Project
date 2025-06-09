#' Differentially Private Proportion Estimator
#'
#' Computes a private estimate of a proportion using uniform random variables.
#'
#' @param theta True proportion
#' @param eps Privacy parameter (ε)
#' @param n Sample size
#'
#' @return Differentially private proportion estimate
get_pi <- function(theta, eps, n) {
  u <- runif(n)
  w <- runif(1, -0.5, 0.5)
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}


#' FIMA algorithm
#'
#' Generates fiducial samples for differentially private inference.
#'
#' @param pi0 Private proportion estimate
#' @param eps Privacy parameter (ε)
#' @param n Sample size
#' @param H Number of fiducial samples (default = 1)
#' @param seed Random seed (default = 123)
#'
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


#' Two-Sample FIMA Test
#'
#' Performs differentially private two-sample proportion comparison.
#'
#' @param pi1 Private proportion estimate for group 1
#' @param pi2 Private proportion estimate for group 2
#' @param H Number of fiducial samples
#' @param eps Privacy parameter (ε)
#' @param n1 Sample size for group 1
#' @param n2 Sample size for group 2
#' @param seed Random seed
#'
#' @return Vector of H differences between group distributions
fima2samples <- function(pi1, pi2, H, eps, n1, n2, seed) {
  set.seed(seed)
  
  distri_grp1 <- fima(pi0 = pi1, H=H, eps = eps/2, n = n1, seed=sample(1:10000,1))
  distri_grp2 <- fima(pi0 = pi2, H=H, eps = eps/2, n = n2, seed=sample(10001:20000,1))
  
  res <- distri_grp1 - distri_grp2
  return(res)
}



#' Simulate Power and Computation Time using FIMA algorithm
#'
#' Evaluates statistical power and computation time across different proportion values.
#'
#' @param H Number of fiducial samples
#' @param B Number of Monte Carlo replications
#' @param n1 Sample size for group 1
#' @param n2 Sample size for group 2
#' @param theta1 True proportion for group 1
#' @param theta2 Vector of true proportions for group 2
#' @param eps Privacy parameter (ε)
#' @param alpha Significance level
#'
#' @return Data frame with:
#'   - theta2: Proportion values for group 2
#'   - power: Estimated power
#'   - mean_time_ms: Mean computation time (ms)
#'   - median_time_ms: Median computation time (ms)
simulate_power_and_time <- function(H,B,n1,n2,theta1,theta2,eps,alpha) {
  # Initialize matrices to store results
  res_PJ2 <- matrix(NA, B, length(theta2))
  time_matrix <- matrix(NA, B, length(theta2))
  
  for (j in 1:length(theta2)) {
    for (i in 1:B) {
      set.seed(i)
      start_time <- Sys.time()
      
      # Private computations
      pi1 <- get_pi(theta = theta1, eps = eps/2, n = n1)
      pi2 <- get_pi(theta = theta2[j], eps = eps/2, n = n2)
      
      # Compute empirical density
      emp_density_PJ <- fima2samples(
        pi1 = pi1, pi2 = pi2, H = H, eps = eps, 
        n1 = n1, n2 = n2, seed = 10 + i
      )
      
      # Store p-value and time
      res_PJ2[i, j] <- (sum(emp_density_PJ >= 0)) / (H + 1)
      time_matrix[i, j] <- as.numeric(Sys.time() - start_time) * 1000  # in ms
      
      if (i %% 1000 == 0) cat(sprintf("Theta2 = %.1f, Replication %d/%d\n", theta2[j], i, B))
    }
    cat(sprintf("Completed Theta2 = %.1f\n", theta2[j]))
  }
  
  # Compute summary statistics
  power <- colMeans(res_PJ2 < alpha, na.rm = TRUE)
  mean_time <- colMeans(time_matrix, na.rm = TRUE)
  median_time <- apply(time_matrix, 2, median, na.rm = TRUE)
  
  # Return as a dataframe
  data.frame(
    theta2 = theta2,
    power = power,
    mean_time_ms = mean_time,
    median_time_ms = median_time
  )
}



# Simulation parameters
H <- 10^4        # Number of fiducial samples
B <- 10^4        # Number of Monte Carlo replications
n1 <- 30         # Group 1 sample size
n2 <- 100        # Group 2 sample size
theta1 <- 0.2    # Fixed proportion for group 1
theta2 <- seq(0.2, 0.9, 0.1)  # Varying proportions for group 2
eps <- 1         # Privacy parameter (ε = 1)
alpha <- 0.05    # Significance level

# Run simulation
power_time_two_sample_over_theta_fima <- simulate_power_and_time(
  H, B, n1, n2, theta1, theta2, eps, alpha
)

# View results
print(power_time_two_sample_over_theta_fima)

# # Plot power curve
# plot(power_time_two_sample_over_theta_fima$theta2, 
#      power_time_two_sample_over_theta_fima$power,
#      type = "b", xlab = "Group 2 Proportion (θ₂)", ylab = "Power",
#      main = "Power Analysis for FIMA Two-Sample Test (θ₁ = 0.2)")

# Save results
#write.csv(power_time_two_sample_over_theta_fima, "fima_power_results.csv", row.names = FALSE)
