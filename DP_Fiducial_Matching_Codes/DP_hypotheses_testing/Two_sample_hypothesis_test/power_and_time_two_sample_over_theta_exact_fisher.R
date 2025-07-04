#' Perform Two-Sample Exact Binomial Test (Fisher's Test)
#'
#' Conducts Fisher's exact test on simulated binomial data and returns both the 
#' p-value and computation time.
#'
#' @param n Sample size for each group
#' @param theta1 True success probability for group 1
#' @param theta2 True success probability for group 2
#'
#' @return A list containing:
#'   - p_value: The p-value from Fisher's exact test (one-sided less)
#'   - time: Computation time in seconds
#'
#' @details
#' 1. Generates binomial samples for two groups
#' 2. Creates 2x2 contingency table
#' 3. Handles edge cases where test would be undefined
#' 4. Performs Fisher's exact test with one-sided alternative
#' 5. Measures and returns computation time
exact_binom_test <- function(n, theta1, theta2) {
  # Generate samples
  x1 <- rbinom(1, size = n, prob = theta1)
  x2 <- rbinom(1, size = n, prob = theta2)
  
  # Create contingency table
  cont_table <- matrix(c(x1, n - x1, x2, n - x2), 
                       nrow = 2,
                       dimnames = list(Group = c("Group1", "Group2"),
                                       Response = c("Success", "Failure")))
  
  # Handle edge cases where test would be undefined
  if ((x1 == 0 && x2 == 0) || (x1 == n && x2 == n)) {
    # These cases perfectly satisfy H0, so p-value should be 1
    return(list(p_value = 1, time = 0))
  }
  
  # Measure computation time
  start_time <- Sys.time()
  # Perform Fisher's exact test (one-sided)
  test_result <- fisher.test(cont_table, alternative = "less")
  end_time <- Sys.time()
  
  return(list(p_value = test_result$p.value, 
              time = as.numeric(end_time - start_time)))
}


#' Simulate Fisher's Exact Test Performance
#'
#' Evaluates the empirical power and computation time of Fisher's 
#' exact test across different true proportion values.
#'
#' @param n Sample size for each group
#' @param theta2_range Vector of true proportion values for group 2
#' @param theta1 True proportion for group 1 (fixed)
#' @param num_simulations Number of Monte Carlo simulations per theta value
#' @param alpha Significance level for empirical power calculation
#'
#' @return A data frame with columns:
#'   - theta: True proportion values for group 2
#'   - empirical_power: Estimated power
#'   - mean_time: Average computation time in milliseconds
#'   - median_time: Median computation time in milliseconds
#'
#' @details
#' For each theta in theta2_range:
#' 1. Runs num_simulations tests comparing theta1 vs theta
#' 2. Calculates empirical power
#' 3. Records computation time statistics
simulate_test_exact <- function(n, theta2_range, theta1, num_simulations, alpha) {
  # Initialize results matrix
  results <- matrix(NA, nrow = length(theta2_range), ncol = 4)
  colnames(results) <- c("theta", "empirical_power", "mean_time", "median_time")
  
  for (i in seq_along(theta2_range)) {
    theta <- theta2_range[i]
    # Simulate under H0 (theta1 = theta2 = theta)
    sim_results <- lapply(1:num_simulations, function(j) {
      set.seed(j)  # Set seed based on simulation number
      exact_binom_test(n, theta1 = theta1, theta2 = theta)
    })
    
    # Extract p-values and times
    p_values <- sapply(sim_results, function(x) x$p_value)
    times <- sapply(sim_results, function(x) x$time)
    
    # Calculate empirical power (proportion of rejections at alpha)
    empirical_power <- mean(p_values < alpha, na.rm = TRUE)
    
    # Calculate mean and median computation time (in milliseconds)
    mean_time <- mean(times, na.rm = TRUE) * 1000
    median_time <- median(times, na.rm = TRUE) * 1000
    
    results[i, ] <- c(theta, empirical_power, mean_time, median_time)
  }
  
  return(as.data.frame(results))
}

# Parameters for simulation
n <- 100                 # Fixed sample size for both groups
theta2_values <- seq(0.2, 0.9, by = 0.1)  # Varying true proportion values for group 2
theta1 <- 0.2            # Fixed proportion for group 1
alpha <- 0.05            # Significance level
num_simulations <- 10^4  # Number of simulations per theta value

# Run simulation
sim_results_exact <- simulate_test_exact(
  n = n, 
  theta2_range = theta2_values, 
  theta1 = theta1, 
  alpha = alpha, 
  num_simulations = num_simulations
)

# View results
print(sim_results_exact)


# Save results
#write.csv(sim_results_exact, "fisher_exact_test_results.csv", row.names = FALSE)
