#' Perform Timed Two-Sample Proportion Test
#'
#' Conducts a two-sample proportion test on simulated binomial data and returns both the 
#' p-value and computation time, with robust handling of edge cases.
#'
#' @param n Sample size for each group
#' @param theta1 True success probability for group 1
#' @param theta2 True success probability for group 2
#'
#' @return A list containing:
#'   - p_value: The p-value from the proportion test (one-sided less)
#'   - time: Computation time in milliseconds
#'
#' @details
#' 1. Generates binomial samples for two groups
#' 2. Handles edge cases where test would be undefined or fail
#' 3. Performs proportion test with one-sided alternative and no continuity correction
#' 4. Measures and returns computation time
prop_test_timed <- function(n, theta1, theta2) {
  # Generate samples
  x1 <- rbinom(1, size = n, prob = theta1)
  x2 <- rbinom(1, size = n, prob = theta2)
  
  # Handle edge cases that would lead to NA p-value
  if (x1 == 0 && x2 == 0) {
    # Both groups have 0 successes - no difference to detect
    p_value <- 1
  } else if (x1 == n && x2 == n) {
    # Both groups have 100% success - no difference to detect
    p_value <- 1
  } else {
    # Time the test execution
    start_time <- Sys.time()
    test_result <- tryCatch({
      prop.test(c(x1, x2), n = c(n, n), 
                alternative = "less", 
                correct = FALSE)
    }, error = function(e) {
      # If prop.test fails, return conservative p-value
      list(p.value = 1)
    })
    end_time <- Sys.time()
    p_value <- ifelse(is.na(test_result$p.value), 1, test_result$p.value)
  }
  
  # Calculate time (only for the actual test cases)
  if (!exists("end_time")) {
    end_time <- start_time <- Sys.time()
  }
  
  return(list(p_value = p_value,
              time = as.numeric(end_time - start_time) * 1000)) # Convert to ms
}



#' Simulate Two-Sample Proportion Test Performance
#'
#' Evaluates the power and computation time of a two-sample proportion test across
#' different sample sizes.
#'
#' @param sample_sizes Vector of sample sizes to evaluate
#' @param theta1 True success probability for group 1
#' @param theta2 True success probability for group 2
#' @param B Number of Monte Carlo replications (default = 1000)
#'
#' @return A data frame with columns:
#'   - sample_size: Sample sizes evaluated
#'   - power: Estimated power at each sample size
#'   - mean_time: Average computation time in milliseconds
#'   - median_time: Median computation time in milliseconds
#'
#' @details
#' For each sample size:
#' 1. Runs B simulations of the proportion test
#' 2. Calculates empirical power (proportion of p-values < 0.05)
#' 3. Records computation time statistics
#' 4. Provides progress updates
simulate_power_timed <- function(sample_sizes, theta1, theta2, B = 1000) {
  results <- data.frame(sample_size = numeric(length(sample_sizes)),
                        power = numeric(length(sample_sizes)),
                        mean_time = numeric(length(sample_sizes)),
                        median_time = numeric(length(sample_sizes)))
  
  for (i in seq_along(sample_sizes)) {
    n <- sample_sizes[i]
    # Run simulations with seed depending on j
    sim_results <- lapply(1:B, function(j) {
      set.seed(j)  # Seed depends on simulation index j
      prop_test_timed(n, theta1, theta2)
    })
    
    # Extract results
    p_values <- sapply(sim_results, function(x) x$p_value)
    times <- sapply(sim_results, function(x) x$time)
    
    # Calculate metrics
    results$sample_size[i] <- n
    results$power[i] <- mean(p_values < 0.05)
    results$mean_time[i] <- mean(times)
    results$median_time[i] <- median(times)
    
    # Print progress
    cat(sprintf("Completed n = %d (%.0f%%)\n", n, i/length(sample_sizes)*100))
  }
  
  return(results)
}



# Simulation parameters
sample_sizes <- c(16, 30, 50, 100, 150, 200, 350, 400, 500) 
B <- 10^4      # Number of Monte Carlo replications
theta1 <- 0.8  # True proportion for group 1
theta2 <- 0.9  # True proportion for group 2

# Run simulation
results <- simulate_power_timed(sample_sizes, theta1, theta2, B)

# View results
print(results)


# Save results
#write.csv(results, "proportion_test_results.csv", row.names = FALSE)
