library(glmnet)
library(nloptr)
library(alabama)
library(parallel)
library(pbmcapply)  

#####Helping Functions######
# Function to compute chi2 stat. based on counts in contingency table
compute_chi2 <- function(contingency_tab){
  n <- sum(contingency_tab)
  marginal_row <- apply(contingency_tab,1,sum)/n 
  maringal_col <- apply(contingency_tab,2,sum)/n
  expected <- outer(marginal_row,maringal_col,FUN="*")*n
  return(sum((contingency_tab-expected)^2/expected))
}
# Laplace noise
laplace_noise <- function(sensitivity, epsilon){
  W <- runif(1, -0.5, 0.5)
  noise <- -sensitivity/epsilon*sign(W)*log(1 - 2*abs(W))
  return(noise)
}
# Laplace mechanism for differential privacy
add_laplace_noise <- function(counts, epsilon) {
  sensitivity <- 2
  noisy_counts <- sapply(counts, function(count) count+laplace_noise(sensitivity, epsilon=epsilon))
  noisy_counts <- pmax(noisy_counts, 0)
  return(haldane_anscombe_correction(noisy_counts))  # Ensure non-negative counts
}
# Haldane anscombe correction for when cells are equal to 0
haldane_anscombe_correction <- function(counts){
  counts[counts==0] <- counts[counts==0]+0.5
  return(counts)
}




################ Gabroadi Method     ##########################





# Laplace random generator
rlaplace <- function(n, scale) {
  u <- runif(n, -0.5, 0.5)
  return(-scale * sign(u) * log(1 - 2 * abs(u)))
}


# 2MLE Function for Algorithm 5: Two Step MLE Calculation
two_step_MLE <- function(w, n, noise_type = "Laplace", gamma = NULL) {
  
  # Step 1: Find the most likely contingency table x_hat given noisy counts w
  
  # Objective function to minimize: either ||w - x||_1 or ||w - x||_2^2 based on noise type
  if (noise_type == "Gaussian") {
    objective <- function(x) {
      return(sum((w - x)^2))
    }
  } else if (noise_type == "Laplace") {
    if (is.null(gamma)) gamma <- 0.01  # A small gamma for Laplace noise
    objective <- function(x) {
      return((1 - gamma) * sum(abs(w - x)) + gamma * sum((w - x)^2))
    }
  } else {
    stop("Unknown noise type! Please use either 'Gaussian' or 'Laplace'.")
  }
  
  
  r <- nrow(w)
  c <- ncol(w)
  
  # Initial guess: distribute n uniformly across the contingency table
  initial_guess <- rep(n / (r * c), r * c)
  
  # Lower bounds for the optimization (non-negativity constraint: x >= 0)
  lower_bounds <- rep(0, r * c)
  
  # Define inequality constraint function: x >= 0 (non-negativity)
  hin <- function(x) {
    return(x)
  }
  
  # Define equality constraint function: sum(x) = n
  heq <- function(x) {
    return(sum(x) - n)
  }
  
  # Run the optimization using constrOptim.nl
  fit <- constrOptim.nl(
    par = initial_guess,    # Initial guess (ensure sum(initial_guess) = n)
    fn = objective,         # Objective function (L1 or L2 norm)
    hin = hin,              # Inequality constraint function (non-negativity)
    heq = heq,              # Equality constraint function (sum(x) = n)
    control.optim = list(maxit = 1000)  # Optional control parameters
  )
  
  # Reshape the optimized solution into the contingency table form
  x_hat <- matrix(fit$par, nrow = r, ncol = c)
  
  # Step 2: Compute the MLE for π^(1) and π^(2) given the denoised contingency table
  
  # MLE for row probabilities π^(1)
  pi_1_hat <- rowSums(x_hat) / n
  
  # MLE for column probabilities π^(2)
  pi_2_hat <- colSums(x_hat) / n
  
  # Rule of thumb: return NULL if any count is less than 5
  if (any(x_hat < 5)) {
    #print(paste("Optimization yields a contingency table with at least one cell having entry that is less than 5", x_hat))
    return(list(pi_1_hat = NULL, pi_2_hat = NULL))
  }
  
  
  # Return the MLE estimates
  return(list(pi_1_hat = pi_1_hat, pi_2_hat = pi_2_hat, x_hat = x_hat))
}



# MC Independence Testing


# Helper function to generate noisy counts (Laplace or Gaussian noise)
generate_noisy_table <- function(x, noise_type = "Laplace", epsilon, delta) {
  r <- nrow(x)
  c <- ncol(x)
  #w <- matrix(0, nrow = r, ncol = c)
  
  if (noise_type == "Laplace") {
    sensitivity = 2
    scale <- sensitivity / epsilon
    noise_matrix <- matrix(rlaplace(r * c, scale = scale), nrow = r, ncol = c)
  } else if (noise_type == "Gaussian") {
    #sensitivity = sqrt(2)
    sigma <- 2 * sqrt(log(2 / delta)) / epsilon
    noise_matrix <- matrix(rnorm(r * c, mean = 0, sd = sigma), nrow = r, ncol = c)
  }
  
  w <- x + noise_matrix
  return(w)
}

# Function to compute chi-squared statistic
compute_chi_squared <- function(x, pi_1_hat, pi_2_hat) {
  n=sum(x)
  expected <- outer(pi_1_hat, pi_2_hat) * n
  chi_sq <- sum((x - expected)^2 / expected)
  return(chi_sq)
}

# MCIndepD algorithm
mc_indep_test <- function(x, epsilon, delta, alpha, noise_type = "Laplace", k = 1000) {
  n <- sum(x)  # Total number of samples
  
  # Step 1: Add noise to the contingency table
  w <- generate_noisy_table(x, noise_type = noise_type, epsilon = epsilon, delta = delta)
  
  # Step 2: Estimate pi_tilde^(1) and pi_tilde^(2) using 2MLE
  mle_result <- two_step_MLE(w, n, noise_type = noise_type)
  pi_1_hat <- mle_result$pi_1_hat
  pi_2_hat <- mle_result$pi_2_hat
  
  
  if (is.null(pi_1_hat) || is.null(pi_2_hat)) {
    return("Fail to Reject H0: The variables are independent.")
  }
  
  # Step 3: Compute the test statistic for observed data
  q_obs <- compute_chi_squared(w, pi_1_hat, pi_2_hat)
  
  # Step 4: Monte Carlo sampling to approximate the threshold
  q_values <- numeric(k)
  for (i in 1:k) {
    # Generate synthetic table from estimated probabilities
    synthetic_x <- matrix(rmultinom(1, n, prob = outer(pi_1_hat, pi_2_hat)), 
                          nrow = length(pi_1_hat), ncol = length(pi_2_hat))
    
    # Add noise to the synthetic table
    w_synthetic <- generate_noisy_table(synthetic_x, noise_type = noise_type, epsilon = epsilon, delta = delta)
    
    # Recompute MLE for the noisy synthetic data
    n_synthetic_x <- sum(synthetic_x)
    mle_synthetic <- two_step_MLE(w_synthetic, n_synthetic_x, noise_type = noise_type)
    if (is.null(mle_synthetic$pi_1) || is.null(mle_synthetic$pi_2)) {
      
      
      next
      
    }
    
    # Compute chi-squared statistic for synthetic data
    q_values[i] <- compute_chi_squared(w_synthetic, mle_synthetic$pi_1_hat, mle_synthetic$pi_2_hat)
  }
  
  
  
  
  # Step 5: Compute the (1 - alpha) quantile of the Monte Carlo statistics
  q_values <- sort(q_values)
  q_threshold <- q_values[ceiling((1 - alpha) * (k+1))]
  
  # Step 6: Make the final decision
  if (q_obs > q_threshold) {
    return("Reject H0: The variables are dependent.")
  } else {
    return("Fail to Reject H0: The variables are independent.")
  }
}



##################  Power/level Analysis ######################




Gabroadi_simulations <- function(joint_probs0, n_vec, epsilon, delta, alpha, 
                                 noise_type = "Laplace", k, sim, num_cores = detectCores()) {
  # Initialize matrices to store results and timings
  results_matrix <- matrix(NA, nrow = length(n_vec), ncol = 2)
  time_matrix <- matrix(NA, nrow = sim, ncol = length(n_vec))
  colnames(results_matrix) <- c("Sample_size", "Gabroadi_power_level")
  
  # Get p_A from joint_probs0 dimensions (assuming square matrix)
  p_A <- sqrt(length(joint_probs0))
  
  for (i in seq_along(n_vec)) {
    n <- n_vec[i]
    
    # Show progress message for current sample size
    message("\nRunning simulations for sample size n = ", n, " (", i, "/", length(n_vec), ")")
    
    # Run simulations in parallel with progress bar
    sim_results <- pbmclapply(1:sim, function(b) {
      set.seed(b * 10 + b)  # Set seed for reproducibility within each core
      start_time <- Sys.time()
      
      result <- tryCatch({
        counts <- rmultinom(1, n, prob = as.vector(joint_probs0))
        contingency_table <- matrix(counts, nrow = p_A, byrow = FALSE)
        
        # Run the Monte Carlo independence test
        test_result <- mc_indep_test(contingency_table, epsilon, delta, alpha, noise_type, k)
        
        # Return both the test result and computation time
        list(
          reject = ifelse(test_result == "Reject H0: The variables are dependent.", 1, 0),
          time = as.numeric(Sys.time() - start_time) * 1000
        )
      }, error = function(e) {
        message("Error in iteration ", b, " for n=", n, ": ", e$message)
        return(list(reject = NA, time = NA))
      })
      
      return(result)
    }, mc.cores = num_cores)
    
    # Extract results from parallel simulations
    reject_counts <- sapply(sim_results, function(x) x$reject)
    time_matrix[, i] <- sapply(sim_results, function(x) x$time)
    
    # Calculate power/level
    valid_results <- reject_counts[!is.na(reject_counts)]
    if (length(valid_results) > 0) {
      results_matrix[i, ] <- c(n, mean(valid_results))
    } else {
      results_matrix[i, ] <- c(n, NA)
    }
  }
  
  # Calculate timing statistics
  mean_times <- apply(time_matrix, 2, mean, na.rm = TRUE)
  median_times <- apply(time_matrix, 2, median, na.rm = TRUE)
  
  # Create results dataframe
  results_df <- data.frame(
    Sample_size = n_vec,
    Power_level = results_matrix[, "Gabroadi_power_level"],
    Mean_time_ms = mean_times,
    Median_time_ms = median_times
  )
  
  return(results_df)
}



##### Type I error simulations example ######
p_A <- c(0.5, 0.5)
p_B <- c(0.5, 0.5)
delta = 0.01
joint_probs0 <-c(0.25, 0.25,0.25, 0.25)  
n_vec <- c(20,30,50, 100, 200, 500, seq(1000, 20000, by = 1000))  
epsilon <- 0.1
k <- 1000
sim <- 10^3
alpha <- 0.05


# Run the simulation with Laplace noise
level_and_time_chisquare_over_n_Gabroadi <-Gabroadi_simulations(joint_probs0 , n_vec, epsilon = epsilon, delta = delta, 
                                                                alpha = alpha, noise_type = "Laplace", k = k, sim = sim)
level_and_time_chisquare_over_n_Gabroadi


### Power #######
joint_probs0 <-c(0.25, 0.25,0.25, 0.25) + delta*c(1,-1,-1,1) 
power_and_time_chisquare_over_n_Gabroadi <- Gabroadi_simulations(joint_probs0 , n_vec, epsilon = epsilon, delta = delta, 
                                                                 alpha = alpha, noise_type = "Laplace", k = k, sim = sim)
power_and_time_chisquare_over_n_Gabroadi
