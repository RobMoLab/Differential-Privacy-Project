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



#' fima
#'
#' Generates fiducial samples for differentially private inference using
#' a beta distribution approach.
#'
#' @param pi0 Private proportion estimate
#' @param eps Privacy parameter (ε)
#' @param n Sample size
#' @param H Number of fiducial samples (default = 1)
#' @param delta Sensitivity parameter (default = 1)
#' @param seed Random seed (default = 123)
#'
#' @return Vector of H fiducial samples
#' @details
#' 1. Adds Laplace noise scaled by privacy budget ε
#' 2. Generates beta-distributed fiducial samples
#' 3. Handles boundary cases with machine epsilon
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





#' Two-Sample Fima Test
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
fima2samples<- function(pi1, pi2, H, eps, n1, n2, seed) {
  set.seed(seed)
  
  distri_grp1 <- fima(pi0 = pi1, H=H, eps = eps/2, n = n1, seed=sample(1:10000,1))
  distri_grp2 <- fima(pi0 = pi2, H=H, eps = eps/2, n = n2, seed=sample(10001:20000,1))
  
  res <- distri_grp1 - distri_grp2
  return(res)
}



#' Compute Power for Two-Sample test with Fima
#'
#' Evaluates statistical power across different sample sizes.
#'
#' @param theta1 True proportion for group 1
#' @param theta2 True proportion for group 2
#' @param n_values Vector of sample sizes to evaluate
#' @param eps Privacy parameter (ε)
#' @param B Number of Monte Carlo replications
#' @param H Number of fiducial samples
#'
#' @return Data frame with:
#'   - n: Sample sizes
#'   - power: Estimated power
#'   - mean_time: Mean computation time (ms)
#'   - median_time: Median computation time (ms)
compute_two_sample_power <- function(theta1, theta2, n_values, eps, B, H) {
  results <- data.frame(
    n = integer(),
    power = numeric(),
    mean_time = numeric(),
    median_time = numeric(),
    stringsAsFactors = FALSE
  )
  
  res_pj2 <- matrix(NA, B, length(n_values))
  compute_times <- matrix(NA, B, length(n_values))
  
  for (j in 1:length(n_values)) {
    n <- n_values[j]
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    
    for (i in 1:B) {
      set.seed(i)
      start_time <- Sys.time()
      
      pi1 <- get_pi(theta = theta1, eps = eps/2, n = n)
      pi2 <- get_pi(theta = theta2, eps = eps/2, n = n)
      
      emp_density_pj <- fima2samples(
        pi1 = pi1, pi2 = pi2, 
        H = H, eps = eps, 
        n1 = n, n2 = n, 
        seed = i + 2*B
      )
      res_pj2[i,j] <- sum(emp_density_pj >= 0)/(H + 1)
      end_time <- Sys.time()
      compute_times[i,j] <- as.numeric(end_time - start_time)*1000
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    results <- rbind(results, data.frame(
      n = n,
      power = mean(res_pj2[,j] < 0.05, na.rm = TRUE),
      mean_time = mean(compute_times[,j], na.rm = TRUE),
      median_time = median(compute_times[,j], na.rm = TRUE)
    ))
    
    print(paste("Completed n =", n))
  }
  
  return(results)
}



# Simulation parameters
H <- 10^4        # Number of fiducial samples
B <- 10^4        # Number of Monte Carlo replications
n_values <- c(16, 30, 50, 100, 150, 200, 350, 400, 500)
theta1 <- 0.8    # Group 1 proportion
theta2 <- 0.9    # Group 2 proportion
eps <- 1         # Privacy parameter (ε = 1)

# Run simulation
results <- compute_two_sample_power(theta1, theta2, n_values, eps, B, H)

# View results
print(results)


# Save results
#write.csv(results, "jini_power_results.csv", row.names = FALSE)
