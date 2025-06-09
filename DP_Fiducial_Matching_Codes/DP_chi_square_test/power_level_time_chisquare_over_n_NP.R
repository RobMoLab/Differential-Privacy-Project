Nonprivate_simulation <- function(joint_probs0, n_vec, alpha, sim) {
  # Initialize matrices to store results and timings
  results_matrix <- matrix(NA, nrow = length(n_vec), ncol = 2)
  time_matrix <- matrix(NA, nrow = sim, ncol = length(n_vec))
  colnames(results_matrix) <- c("Sample_Size", "Power")
  
  for (i in seq_along(n_vec)) {
    n <- n_vec[i]
    p_values <- numeric(sim)  # Store p-values for current sample size
    
    for (b in 1:sim) {
      set.seed(b * 10 + b)
      start_time <- Sys.time()
      
      # Generate data and perform test
      counts <- rmultinom(1, n, prob = joint_probs0)
      nonprivate_contingency_table <- matrix(counts, nrow = 2)
      test_result <- tryCatch({
        chisq.test(nonprivate_contingency_table)
      }, error = function(e) {
        message("Error in iteration ", b, " for n=", n, ": ", e$message)
        return(list(p.value = NA))
      })
      
      # Store results and timing
      p_values[b] <- test_result$p.value
      time_matrix[b, i] <- as.numeric(Sys.time() - start_time) * 1000  # ms
    }
    
    # Calculate power for current sample size
    NonPrivate_power <- mean(p_values <= alpha, na.rm = TRUE)
    results_matrix[i, ] <- c(n, NonPrivate_power)
  }
  
  # Calculate timing statistics (handle NAs)
  mean_times <- apply(time_matrix, 2, mean, na.rm = TRUE)
  median_times <- apply(time_matrix, 2, median, na.rm = TRUE)
  
  # Replace any remaining NA/NaN with 0 for plotting
  mean_times[is.na(mean_times)] <- 0
  median_times[is.na(median_times)] <- 0
  
  # Create results dataframe
  results_df <- data.frame(
    Sample_Size = n_vec,
    Power = results_matrix[, "Power"],
    Mean_Time_ms = mean_times,
    Median_Time_ms = median_times
  )
  
  return(list(results = results_df, time_matrix = time_matrix))
}

# Level
delta = 0.01
joint_probs0 <-c(0.25, 0.25,0.25, 0.25) 
n_vec <- c(20,30,50, 100, 200, 500, seq(1000, 20000, by = 1000))
sim <- 10^3 
alpha <- 0.05

level_time_chisquare_over_n_NP <- Nonprivate_simulation(joint_probs0, n_vec, alpha, sim)

print(level_time_chisquare_over_n_NP$results)


# power
joint_probs0 <-c(0.25, 0.25,0.25, 0.25) + delta*c(1,-1,-1,1) 

power_time_chisquare_over_n_NP <- Nonprivate_simulation(joint_probs0, n_vec, alpha, sim)


print(power_time_chisquare_over_n_NP$results)
