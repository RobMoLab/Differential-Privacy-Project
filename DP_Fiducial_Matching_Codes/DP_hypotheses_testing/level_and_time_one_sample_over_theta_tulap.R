

#' Jordan's functions (FU)
#' 
#' Helper function for Tulap distribution calculations. Maps values from [-½,½] to [0,1].
#' 
#' @param u Numeric input in range [-0.5, 0.5]
#' @return Transformed value in [0,1]
FU <- function(u) {
  return(u + 1/2)
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


#' Compute Type I Error for DP Binomial Test
#' 
#' Evaluates Type I error rates of a differentially private binomial test
#' using Tulap noise across different null proportions.
#' 
#' @param theta0 Vector of null proportions to evaluate
#' @param n Sample size
#' @param eps Privacy parameter (ε)
#' @param B Number of Monte Carlo replications
#' 
#' @return Data frame containing:
#'   - theta: Null proportion values
#'   - type1_error: Estimated Type I error rate
#'   - mean_time: Average computation time (ms)
#'   - median_time: Median computation time (ms)
#' 
#' @details
#' For each θ₀ in theta0:
#' 1. Generates B binomial samples from Binom(n, θ₀)
#' 2. Adds Tulap noise for differential privacy
#' 3. Computes p-values using Tulap CDF
#' 4. Estimates Type I error rate as proportion of p-values < 0.05
#' 5. Tracks computation times
compute_type1_error <- function(theta0, n, eps, B) {
  # Initialize results structure
  results <- data.frame(
    theta = numeric(),
    type1_error = numeric(),
    mean_time = numeric(),
    median_time = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Initialize storage matrices
  res_tulap <- matrix(NA, B, length(theta0))
  compute_times <- matrix(NA, B, length(theta0))
  
  # Main simulation loop
  for (j in 1:length(theta0)) {
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    
    for (i in 1:B) {
      set.seed(i)
      start_time <- Sys.time()
      
      # Generate data under null
      X <- rbinom(n = 1, size = n, prob = theta0[j])
      values <- seq(0, n)
      pdf <- dbinom(values, size = n, prob = theta0[j])
      
      # Generate Tulap noise
      U <- runif(n = 1, min = -1/2, max = 1/2)  # Uniform component
      G1 <- rgeom(n = 1, prob = 1 - exp(-eps))  # Geometric components
      G2 <- rgeom(n = 1, prob = 1 - exp(-eps))
      Tulap <- X + U + G1 - G2  # Noisy statistic
      
      # Compute p-value using Tulap CDF
      cdf <- ptulap(values - Tulap, exp(-eps))
      res_tulap[i,j] <- t(cdf) %*% pdf
      
      # Record computation time (ms)
      compute_times[i,j] <- as.numeric(Sys.time() - start_time) * 1000
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Store results for this theta0
    results <- rbind(results, data.frame(
      theta = theta0[j],
      type1_error = mean(res_tulap[,j] < 0.05),
      mean_time = mean(compute_times[,j]),
      median_time = median(compute_times[,j])
    ))
    
    print(paste("Completed theta =", theta0[j]))
  }
  
  return(results)
}

# Simulation parameters
B <- 10^4       # Number of replications
n <- 100        # Sample size
theta0 <- seq(0.1, 0.9, by = 0.1)  # Null proportions to evaluate
eps <- 1        # Privacy parameter 

# Run simulation
results <- compute_type1_error(theta0, n, eps, B)

# View results
print(results)

# Optional: Save results
# write.csv(results, "dp_type1_error_results.csv", row.names = FALSE)
