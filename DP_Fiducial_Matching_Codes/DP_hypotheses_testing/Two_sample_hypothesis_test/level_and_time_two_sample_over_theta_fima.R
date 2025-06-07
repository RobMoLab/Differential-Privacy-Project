#' Private Proportion Estimator
#'
#' Computes a differentially private estimate of a proportion.
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


#' FIMA
#'
#' Generates fiducial samples for differentially private inference.
#' 
#' @param pi0 Private proportion estimate
#' @param eps Privacy parameter (ε)
#' @param n Sample size
#' @param H Number of fiducial samples (default = 1)
#' @param delta Sensitivity parameter (default = 1)
#' @param seed Random seed (default = 123)
#' 
#' @return Vector of H fiducial samples
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
    n*pi_star[index] + 0.5,  # Haldane-Anscombe correction
    n*(1 - pi_star[index]) + 0.5
  )
  
  # Handle boundary cases
  theta_ifb[pi_star >= 1] <- 1 - .Machine$double.eps
  theta_ifb[pi_star <= 0] <- .Machine$double.eps
  
  return(theta_ifb)
}


#' Two-Sample FIMA Test
#'
#' Computes the difference between two groups using FIMA.
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
fima2sample <- function(pi1, pi2, H, eps, n1, n2, seed) {
  set.seed(seed)
  
  # Generate distributions for both groups
  distri_grp1 <- fima(pi0 = pi1, H=H, eps = eps/2, n = n1, seed=sample(1:10000,1))
  distri_grp2 <- fima(pi0 = pi2, H=H, eps = eps/2, n = n2, seed=sample(10001:20000,1))
  
  # Compute differences
  res <- distri_grp1 - distri_grp2
  return(res)
}




#' Compute Type I Error for Two-Sample FIMA Test
#'
#' Evaluates Type I error rates under the null hypothesis θ₁ = θ₂.
#' 
#' @param theta0 Vector of null proportion values
#' @param n1 Sample size for group 1
#' @param n2 Sample size for group 2
#' @param eps Privacy parameter (ε)
#' @param B Number of Monte Carlo replications
#' @param H Number of fiducial samples
#' @param alpha Significance level (default = 0.05)
#' 
#' @return Data frame with:
#'   - theta: Proportion values
#'   - type1_error: Estimated Type I error rates
#'   - mean_time: Mean computation time (ms)
#'   - median_time: Median computation time (ms)
compute_type1_error <- function(theta0, n1, n2, eps, B, H, alpha = 0.05) {
  results <- data.frame(
    theta = numeric(),
    type1_error = numeric(),
    mean_time = numeric(),
    median_time = numeric(),
    stringsAsFactors = FALSE
  )
  
  res_pj2 <- matrix(NA, B, length(theta0))
  compute_times <- matrix(NA, B, length(theta0))
  
  for (j in 1:length(theta0)) {
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    
    for (i in 1:B) {
      set.seed(i)
      start_time <- Sys.time()
      
      # Generate private estimates
      pi1 <- get_pi(theta = theta0[j], eps = eps/2, n = n1)
      pi2 <- get_pi(theta = theta0[j], eps = eps/2, n = n2)
      
      # Compute p-value
      emp_density_pj <- fima2sample(
        pi1 = pi1, pi2 = pi2, 
        H = H, eps = eps, 
        n1 = n1, n2 = n2, 
        seed = i + 2*B
      )
      res_pj2[i,j] <- sum(emp_density_pj < 0)/(H + 1)
      
      compute_times[i,j] <- as.numeric(Sys.time() - start_time)*1000
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    results <- rbind(results, data.frame(
      theta = theta0[j],
      type1_error = mean(res_pj2[,j] < alpha, na.rm = TRUE),
      mean_time = mean(compute_times[,j], na.rm = TRUE),
      median_time = median(compute_times[,j], na.rm = TRUE)
    ))
    print(paste("Completed theta =", theta0[j]))
  }
  return(results)
}



# Simulation parameters
H <- 10^4    # Number of fiducial samples
B <- 10^4    # Number of Monte Carlo replications
n1 <- 30     # Group 1 sample size
n2 <- 30     # Group 2 sample size
theta0 <- seq(0.1, 0.9, by = 0.1)  # Null proportions to evaluate
eps <- 1     # Privacy parameter (ε = 1)
alpha <- 0.05 # Significance level

# Run simulation
sim_results <- compute_type1_error(theta0, n1, n2, eps, B, H, alpha)

# View results
print(sim_results)


# Save results
#write.csv(sim_results, "fima_two_sample_results.csv", row.names = FALSE)
