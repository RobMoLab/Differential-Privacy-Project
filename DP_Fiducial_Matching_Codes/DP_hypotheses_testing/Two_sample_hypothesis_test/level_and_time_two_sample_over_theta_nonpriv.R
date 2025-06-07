#' Perform Two-Sample Proportion Test
#'
#' Conducts a two-sample proportion test on simulated binomial data and returns both the 
#' p-value and computation time. Uses prop.test() without continuity correction.
#'
#' @param n Sample size for each group
#' @param theta1 True success probability for group 1
#' @param theta2 True success probability for group 2
#'
#' @return A list containing:
#'   - p_value: The p-value from the proportion test (one-sided less)
#'   - time: Computation time in seconds
#'
#' @details
#' 1. Generates binomial samples for two groups
#' 2. Handles edge cases where test would be undefined (all successes/failures)
#' 3. Performs proportion test with one-sided alternative and no continuity correction
#' 4. Measures and returns computation time
prop_test <- function(n, theta1, theta2) {
  # Generate samples
  x1 <- rbinom(1, size = n, prob = theta1)
  x2 <- rbinom(1, size = n, prob = theta2)
  
  # Handle edge cases where test would be undefined
  if ((x1 == 0 && x2 == 0) || (x1 == n && x2 == n)) {
    # These cases perfectly satisfy H0, so p-value should be 1
    return(list(p_value = 1, time = 0))
  }
  
  # Measure computation time
  start_time <- Sys.time()
  # Perform two-sample proportion test (one-sided)
  test_result <- prop.test(c(x1, x2), n = c(n, n), 
                           alternative = "less", 
                           correct = FALSE)
  end_time <- Sys.time()
  
  return(list(p_value = test_result$p.value, 
              time = as.numeric(end_time - start_time)))
}



#' Simulate Two-Sample Proportion Test Performance
#'
#' Evaluates the empirical Type I error rate and computation time of a two-sample
#' proportion test across different true proportion values.
#'
#' @param n Sample size for each group
#' @param theta_range Vector of true proportion values to evaluate
#' @param num_simulations Number of Monte Carlo simulations per theta value
#' @param alpha Significance level for Type I error calculation (default = 0.05)
#'
#' @return A data frame with columns:
#'   - theta: True proportion value
#'   - empirical_level: Estimated Type I error rate
#'   - mean_time: Average computation time in milliseconds
#'   - median_time: Median computation time in milliseconds
#'
#' @details
#' For each theta in theta_range:
#' 1. Runs num_simulations tests under H0 (theta1 = theta2 = theta)
#' 2. Calculates empirical Type I error rate
#' 3. Records computation time statistics
simulate_test <- function(n, theta_range, num_simulations, alpha=0.05) {
  # Initialize results matrix
  results <- matrix(NA, nrow = length(theta_range), ncol = 4)
  colnames(results) <- c("theta", "empirical_level", "mean_time", "median_time")
  
  for (i in seq_along(theta_range)) {
    theta <- theta_range[i]
    # Simulate under H0 (theta1 = theta2 = theta)
    sim_results <- lapply(1:num_simulations, function(j) {
      set.seed(j)  # Set seed based on simulation number
      prop_test(n, theta1 = theta, theta2 = theta)
    })
    
    # Extract p-values and times
    p_values <- sapply(sim_results, function(x) x$p_value)
    times <- sapply(sim_results, function(x) x$time)
    
    # Calculate empirical level (proportion of rejections at alpha)
    empirical_level <- mean(p_values < alpha, na.rm = TRUE)
    
    # Calculate mean and median computation time (in milliseconds)
    mean_time <- mean(times, na.rm = TRUE) * 1000
    median_time <- median(times, na.rm = TRUE) * 1000
    
    results[i, ] <- c(theta, empirical_level, mean_time, median_time)
  }
  
  return(as.data.frame(results))
}





# Parameters for simulation
sample_size <- 30          # Fixed sample size for both groups
theta_values <- seq(0.1, 0.9, by = 0.1)  # True proportion values to evaluate
alpha <- 0.05             # Significance level
num_simulations <- 10^4    # Number of simulations per theta value

# Run simulation
sim_results <- simulate_test(
  n = sample_size, 
  theta_range = theta_values, 
  num_simulations = num_simulations
)

# View results
print(sim_results)


# Save results
#write.csv(sim_results, "proportion_test_results.csv", row.names = FALSE)
