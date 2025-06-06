#' Private Proportion Estimator
#' 
#' Computes a differentially private estimate of a proportion using Laplace mechanism.
#' 
#' @param theta True proportion
#' @param u Vector of uniform random variables (length n)
#' @param w Uniform random variable in [-0.5, 0.5] for noise generation
#' @param eps Privacy parameter (ε)
#' @param n Sample size
#' 
#' @return Differentially private estimate of the proportion
get_pi <- function(theta, u, w, eps, n) {
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}


#' FIMA Algorithm
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
#' @return Vector of fiducial samples
fima <- function(pi0, eps, n, H = 1, delta = 1, seed = 123) {
  set.seed(seed)
  
  # Generate noise variables
  Wj <- runif(H, -0.5, 0.5)
  
  # Create perturbed proportions
  pi_star <- pi0 + delta/(eps*n)*sign(Wj)*log(1 - 2*abs(Wj))
  
  # Generate fiducial samples
  theta_ifb <- rep(NA, H)
  
  # Handle different cases for pi_star values
  index <- pi_star < 1 & pi_star > 0
  theta_ifb[index] <- rbeta(
    sum(index), 
    n*pi_star[index] + 0.5, 
    n*(1 - pi_star[index]) + 0.5
  )
  
  # Boundary cases
  theta_ifb[pi_star >= 1] <- 1 - .Machine$double.eps
  theta_ifb[pi_star <= 0] <- .Machine$double.eps
  
  return(theta_ifb)
}


#' Power Simulation for FIMA Test
#' 
#' Evaluates the power of the FIMA-based hypothesis test across sample sizes.
#' 
#' @param theta0 True proportion under alternative
#' @param n_values Vector of sample sizes to evaluate
#' @param eps Privacy parameter (ε)
#' @param B Number of Monte Carlo replications
#' @param H Number of fiducial samples
#' 
#' @return Data frame with power and timing results
compute_power <- function(theta0, n_values, eps, B, H) {
  # Initialize results data structure
  results <- data.frame(
    n = integer(),
    power_fima = numeric(),
    mean_time_fima = numeric(),
    median_time_fima = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Initialize storage matrices
  res_jini <- matrix(NA, B, length(n_values))
  compute_times <- matrix(NA, B, length(n_values))
  
  # Main simulation loop
  for (j in 1:length(n_values)) {
    n <- n_values[j]
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    
    for (i in 1:B) {
      set.seed(i)
      start_time <- Sys.time()
      
      # Generate private estimate
      pi0 <- get_pi(
        theta = theta0,
        u = runif(n), 
        w = runif(1, -0.5, 0.5),
        eps = eps,
        n = n
      )
      
      # Compute FIMA p-value
      emp_density <- fima(
        pi0 = pi0,
        H = H,
        eps = eps,
        n = n,
        seed = i + 2*B
      )
      res_jini[i,j] <- sum(emp_density <= 0.9)/(H + 1)
      
      # Record computation time (ms)
      end_time <- Sys.time()
      compute_times[i,j] <- as.numeric(end_time - start_time)*1000
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Store results for this sample size
    results <- rbind(results, data.frame(
      n = n,
      power_fima = mean(res_jini[,j] < 0.05),
      mean_time_fima = mean(compute_times[,j]),
      median_time_fima = median(compute_times[,j])
    ))
    
    print(paste("Completed n =", n))
  }
  
  return(results)
}


# Simulation parameters
H <- 10^3    # Number of fiducial samples
B <- 10^4    # Number of Monte Carlo replications
n_values <- c(16, 30, 50, 100, 150, 200, 350, 400, 500)
theta0 <- 0.95  # True proportion under alternative
eps <- 1        # Privacy parameter (ε = 1)

# Run simulation
sim_results <- compute_power(theta0, n_values, eps, B, H)

# View results
print(sim_results)

# Optional: Save results
# write.csv(sim_results, "fima_power_results.csv", row.names = FALSE)
