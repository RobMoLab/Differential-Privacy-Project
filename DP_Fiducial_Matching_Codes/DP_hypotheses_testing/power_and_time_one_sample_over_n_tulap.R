######### Jordan's functions
FU = function(u){# domain is -1/2 to 1/2
  return(u+1/2)
}



#' Tulap Cumulative Distribution Function
#' 
#' Computes the CDF of the Tulap (Truncated-Uniform-Laplace) distribution.
#' 
#' @param t Evaluation points
#' @param b Parameter controlling noise level (b = exp(-ε))
#' @return Vector of CDF values

ptulap = function(t,b){
  cdf = function(t,b){
    ifelse(t<=0,b^(-round(t))/(1+b)*(b+FU(t-round(t))*(1-b)),1-b^(round(t))/(1+b)*(b+FU(round(t)-t)*(1-b)))
  }
  return(sapply(t,cdf,b=b))
}



#' Simulate Power of Differentially Private Binomial Test
#' 
#' Evaluates the power of a Tulap-noise-based DP test for H₀: θ ≤ θ₀ vs Hₐ: θ > θ₀.
#' 
#' @param B Number of Monte Carlo replications
#' @param n_vec Vector of sample sizes to evaluate
#' @param theta0 True proportion under null hypothesis
#' @param eps Privacy parameter (ε)
#' 
#' @return Data frame containing:
#'   - n: Sample sizes
#'   - power_tulap: Estimated power at each sample size
#'   - mean_time_tulap: Average computation time (ms)
simulate_tulap_power <- function(B, n_vec, theta0, eps) {
  # Initialize result matrices
  res_tulap <- matrix(NA, B, length(n_vec))
  time_tulap <- matrix(NA, B, length(n_vec))
  
  # Simulation loop across sample sizes
  for (j in 1:length(n_vec)) {
    for (i in 1:B) {
      set.seed(i)  # Ensure reproducibility
      start_time <- Sys.time()
      
      ###### Generate data and compute Tulap P-value
      # Generate binomial data under alternative (θ = theta0)
      X <- rbinom(n = 1, size = n_vec[j], prob = theta0)
      
      # Compute null distribution (θ = 0.9)
      values <- seq(0, n_vec[j])
      pdf <- dbinom(values, size = n_vec[j], prob = 0.9)
      
      # Generate Tulap noise
      U <- runif(n = 1, min = -1/2, max = 1/2)  # Uniform component
      G1 <- rgeom(n = 1, prob = 1 - exp(-eps))  # Geometric components
      G2 <- rgeom(n = 1, prob = 1 - exp(-eps))
      Tulap <- X + U + G1 - G2  # Noisy statistic
      
      # Compute p-value using Tulap CDF
      cdf <- ptulap(values - Tulap, exp(-eps))
      res_tulap[i,j] <- t(cdf) %*% pdf  # p-value
      
      # Record computation time (ms)
      time_tulap[i,j] <- as.numeric(Sys.time() - start_time) * 1000
    }
  }
  
  # Calculate power (proportion of p-values < 0.05)
  power <- apply(res_tulap < 0.05, 2, mean)
  mean_time <- apply(time_tulap, 2, mean)
  
  # Return results as data frame
  results <- data.frame(
    n = n_vec,
    power_tulap = power,
    mean_time_tulap = mean_time
  )
  
  return(results)
}



# Simulation parameters
B <- 10^4         # Number of Monte Carlo replications
n_vec <- c(16, 30, 50, 100, 150, 200, 350, 400, 500)  # Sample sizes
theta0 <- 0.95    # True proportion under alternative
eps <- 1          # Privacy parameter (ε = 1)

# Run simulation
results <- simulate_tulap_power(B = B, n_vec = n_vec, theta0 = theta0, eps = eps)

# View results
print(results)

# Optional: Save results
# write.csv(results, "tulap_power_results.csv", row.names = FALSE)
