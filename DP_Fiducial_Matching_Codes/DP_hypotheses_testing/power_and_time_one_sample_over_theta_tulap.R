######### Jordan's functions
FU = function(u){# domain is -1/2 to 1/2
  return(u+1/2)
}

ptulap = function(t,b){
  cdf = function(t,b){
    ifelse(t<=0,b^(-round(t))/(1+b)*(b+FU(t-round(t))*(1-b)),1-b^(round(t))/(1+b)*(b+FU(round(t)-t)*(1-b)))
  }
  return(sapply(t,cdf,b=b))
}

#' Compute Empirical Power and Timing for Differentially Private Binomial Test using Tulap method

#' @param theta0 Numeric: Null hypothesis probability (must be in [0, 1]).
#' @param theta_values Numeric vector: Alternative probabilities to test (must be in [0, 1]).
#' @param n Integer: Sample size (must be > 0).
#' @param eps Numeric: Privacy budget (must be > 0).
#' @param B Integer: Number of Monte Carlo replications (must be > 0).
#' @param alpha Numeric: Significance level (default = 0.05, must be in (0, 1)).
#'
#' @return A `data.frame` with columns:
#'   - `theta0`: Null hypothesis value.
#'   - `theta`: Alternative hypothesis value.
#'   - `power_tulap`: Empirical power (proportion of p-values < `alpha`).
#'   - `mean_time_tulap`: Mean computation time (milliseconds).
#'   - `median_time_tulap`: Median computation time (milliseconds).
#'
#' @details
#' For each `theta` in `theta_values`:
#' 1. Simulates `B` binomial samples `X ~ Binom(n, theta)`.
#' 2. Adds Tulap noise: `Tulap = X + U + G1 - G2`.
#' 3. Computes p-values using the Tulap CDF.
#' 4. Estimates power as the proportion of p-values < `alpha`.
#'
#' @examples
#' theta0 <- 0.2
#' theta_values <- seq(0.2, 0.9, by = 0.1)
#' n <- 100
#' eps <- 1
#' B <- 100
#' power_results <- compute_power(theta0, theta_values, n, eps, B)


# Function to compute empirical power with timing
compute_power <- function(theta0, theta_values, n, eps, B, alpha = 0.05) {
  
  # Input validation
  stopifnot(
    theta0 >= 0 && theta0 <= 1,
    all(theta_values >= 0 & theta_values <= 1),
    n > 0,
    eps > 0,
    B > 0,
    alpha > 0 && alpha < 1
  )
  results <- data.frame(
    theta0 = numeric(),
    theta = numeric(),
    power = numeric(),
    mean_time = numeric(),
    median_time = numeric(),
    stringsAsFactors = FALSE
  )
  
  res_tulap <- matrix(NA, B, length(theta_values))
  compute_times <- matrix(NA, B, length(theta_values))
  
  for (j in 1:length(theta_values)) {
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    
    for (i in 1:B) {
      set.seed(i)
      start_time <- Sys.time()
      
      # Generate data from alternative distribution (theta_values[j])
      X <- rbinom(n = 1, size = n, prob = theta_values[j])
      values <- seq(0, n)
      
      # Calculate p-value under H0 (using theta0)
      pdf <- dbinom(values, size = n, prob = theta0)
      
      # Tulap mechanism
      U <- runif(n = 1, min = -1/2, max = 1/2)
      G1 <- rgeom(n = 1, prob = 1 - exp(-eps))
      G2 <- rgeom(n = 1, prob = 1 - exp(-eps))
      Tulap <- X + U + G1 - G2
      cdf <- ptulap(values - Tulap, exp(-eps))
      res_tulap[i,j] <- t(cdf)%*%pdf
      
      compute_times[i,j] <- as.numeric(Sys.time() - start_time)*1000
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Calculate power (probability of rejecting H0 when Ha is true)
    results <- rbind(results, data.frame(
      theta0 = theta0,
      theta = theta_values[j],
      power_tulap = mean(res_tulap[,j] < alpha),
      mean_time_tulap = mean(compute_times[,j]),
      median_time_tulap = median(compute_times[,j])
    ))
    
    print(paste("Completed theta =", theta_values[j]))
  }
  
  return(results)
}




# Simulation parameters
theta0 <- 0.2  # Null hypothesis
theta_values <- seq(0.2, 0.9, by = 0.1)  
n <- 100       # Sample size
eps <- 1       # Privacy budget
B <- 10^4       # Monte Carlo replications
alpha <- 0.05  # Significance level

# Run simulation
power_results <- compute_power(theta0, theta_values, n, eps, B, alpha)

