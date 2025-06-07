#' Perform Timed Fisher's Exact Test
#'
#' Conducts Fisher's exact test on simulated binomial data and returns both the 
#' p-value and computation time, with robust handling of edge cases.
#'
#' @param n Sample size for each group
#' @param theta1 True success probability for group 1
#' @param theta2 True success probability for group 2
#'
#' @return A list containing:
#'   - p_value: The p-value from Fisher's exact test (one-sided less)
#'   - time: Computation time in milliseconds
#'
#' @details
#' 1. Generates binomial samples for two groups
#' 2. Creates 2x2 contingency table
#' 3. Handles edge cases where test would be undefined or fail
#' 4. Performs Fisher's exact test with one-sided alternative
#' 5. Measures and returns computation time
exact_test_timed <- function(n, theta1, theta2) {
  # Generate samples
  x1 <- rbinom(1, size = n, prob = theta1)
  x2 <- rbinom(1, size = n, prob = theta2)
  
  # Create 2x2 contingency table
  cont_table <- matrix(c(x1, n - x1, x2, n - x2), 
                       nrow = 2,
                       dimnames = list(Group = c("Group1", "Group2"),
                                       Response = c("Success", "Failure")))
  
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
      fisher.test(cont_table, alternative = "less")
    }, error = function(e) {
      # If fisher.test fails, return conservative p-value
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


#' Perform Timed Fisher's Exact Test
#'
#' Conducts Fisher's exact test on simulated binomial data and returns both the 
#' p-value and computation time, with robust handling of edge cases.
#'
#' @param n Sample size for each group
#' @param theta1 True success probability for group 1
#' @param theta2 True success probability for group 2
#'
#' @return A list containing:
#'   - p_value: The p-value from Fisher's exact test (one-sided less)
#'   - time: Computation time in milliseconds
#'
#' @details
#' 1. Generates binomial samples for two groups
#' 2. Creates 2x2 contingency table
#' 3. Handles edge cases where test would be undefined or fail
#' 4. Performs Fisher's exact test with one-sided alternative
#' 5. Measures and returns computation time
exact_test_timed <- function(n, theta1, theta2) {
  # Generate samples
  x1 <- rbinom(1, size = n, prob = theta1)
  x2 <- rbinom(1, size = n, prob = theta2)
  
  # Create 2x2 contingency table
  cont_table <- matrix(c(x1, n - x1, x2, n - x2), 
                       nrow = 2,
                       dimnames = list(Group = c("Group1", "Group2"),
                                       Response = c("Success", "Failure")))
  
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
      fisher.test(cont_table, alternative = "less")
    }, error = function(e) {
      # If fisher.test fails, return conservative p-value
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



# Simulation parameters
sample_sizes <- c(16, 30, 50, 100, 150, 200, 350, 400, 500) 
B <- 10^4      # Number of Monte Carlo replications
theta1 <- 0.8  # True proportion for group 1
theta2 <- 0.9  # True proportion for group 2

# Run simulation
results_exact <- simulate_power_timed_exact(sample_sizes, theta1, theta2, B)

# View results
print(results_exact)



# Save results
#write.csv(results_exact, "fisher_exact_test_results.csv", row.names = FALSE)
