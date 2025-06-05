#' Compute Chi-squared Statistic for Contingency Table
#' 
#' @param contingency_tab A 2x2 contingency table of counts
#' @return The chi-squared test statistic
compute_chi2 <- function(contingency_tab){
  n <- sum(contingency_tab)
  marginal_row <- apply(contingency_tab,1,sum)/n 
  maringal_col <- apply(contingency_tab,2,sum)/n
  expected <- outer(marginal_row,maringal_col,FUN="*")*n
  return(sum((contingency_tab-expected)^2/expected))
}

#' Generate Laplace Noise for Differential Privacy
#' 
#' @param sensitivity The sensitivity of the query
#' @param epsilon The privacy budget parameter
#' @return Laplace noise value
laplace_noise <- function(sensitivity, epsilon){
  W <- runif(1, -0.5, 0.5)
  noise <- -sensitivity/epsilon*sign(W)*log(1 - 2*abs(W))
  return(noise)
}

#' Apply Laplace Mechanism to Count Data
#' 
#' @param counts Vector of counts to privatize
#' @param epsilon Privacy budget parameter
#' @return Noisy counts with non-negative values
add_laplace_noise <- function(counts, epsilon) {
  sensitivity <- 2
  noisy_counts <- sapply(counts, function(count) count+laplace_noise(sensitivity, epsilon=epsilon))
  noisy_counts <- pmax(noisy_counts, 0)
  return(haldane_anscombe_correction(noisy_counts))
}

#' Haldane-Anscombe Correction for Zero Counts
#' 
#' Adds 0.5 to zero counts to avoid computational issues
#' @param counts Vector of counts
#' @return Adjusted counts
haldane_anscombe_correction <- function(counts){
  counts[counts==0] <- counts[counts==0]+0.5
  return(counts)
}

#' Privatize Proportions with Differential Privacy
#' 
#' @param prop Original proportion
#' @param n Sample size
#' @param epsilon Privacy budget
#' @return Noisy proportion
priv_prop <- function(prop, n, epsilon){
  W <- runif(1, -0.5, 0.5)
  noisy_prop <-prop -2/(n*epsilon)*sign(W)*log(1 - 2*abs(W))
  return(noisy_prop)
}









#' Simulate Power/Level of Differentially Private Chi-squared Test
#' 
#' @param joint_probs0 True joint probabilities for the null distribution
#' @param n_vec Vector of sample sizes to evaluate
#' @param epsilon Privacy budget
#' @param alpha Significance level
#' @param H Number of bootstrap samples
#' @param sim Number of Monte Carlo simulations
#' @return Data frame with power and computation time results
fima_simulation <- function(joint_probs0, n_vec, epsilon, alpha, H, sim) {
  # Initialize results matrices
  power_results <- matrix(NA, nrow = length(n_vec), ncol = 2)
  time_results <- matrix(NA, nrow = length(n_vec), ncol = 3)
  colnames(power_results) <- c("Sample_size", "Private_power")
  colnames(time_results) <- c("Sample_size", "Mean_time", "Median_time")
  
  for (i in seq_along(n_vec)) {
    n <- n_vec[i]
    
    # Parallel computation across simulations
    results <- mclapply(1:sim, function(b) {
      set.seed(b*10 + b)
      start_time <- Sys.time()
      
      # Generate data and privatize
      counts <- rmultinom(1, n, prob = joint_probs0)
      noisy_contingency_table <- matrix(add_laplace_noise(counts, epsilon), nrow = 2)
      chi2_0b <- compute_chi2(noisy_contingency_table)
      
      # Bootstrap under null
      denoised_counts <- noisy_contingency_table
      n_star <- sum(denoised_counts)
      marginal_row <- apply(denoised_counts, 1, sum) / n_star
      marginal_col <- apply(denoised_counts, 2, sum) / n_star
      
      # Vectorized beta sampling for marginals
      marginal_row_tilde <- matrix(rbeta(H*length(marginal_row), n_star*marginal_row + 0.5, n_star*(1 - marginal_row) + 0.5),
                                   nrow = H, byrow = TRUE)
      marginal_col_tilde <- matrix(rbeta(H*length(marginal_col), n_star*marginal_col + 0.5, n_star*(1 - marginal_col) + 0.5),
                                   nrow = H, byrow = TRUE)
      
      # Generate joint probabilities and noisy counts
      joint_probs_tilde <- array(NA, dim = c(H, 2, 2))
      for (h in 1:H) {
        joint_probs_tilde[h,,] <- outer(marginal_row_tilde[h,], marginal_col_tilde[h,])
      }
      
      counts_tilde_star <- t(sapply(1:H, function(h) {
        rmultinom(1, n_star, prob = joint_probs_tilde[h,,])
      }))
      counts_tilde_star <- t(apply(counts_tilde_star, 1, haldane_anscombe_correction))
      noisy_counts_tilde <- t(apply(counts_tilde_star, 1, function(x) add_laplace_noise(x, epsilon)))
      chi2_hb <- apply(noisy_counts_tilde, 1, function(x) compute_chi2(matrix(x, nrow = 2)))
      
      pval <- mean(chi2_hb >= chi2_0b)
      end_time <- Sys.time()
      
      return(list(pval = pval, time = as.numeric(end_time - start_time) * 1000)) # in milliseconds
    }, mc.cores = detectCores())
    
    # Aggregate results
    p_values <- sapply(results, function(x) x$pval)
    computation_times <- sapply(results, function(x) x$time)
    
    Private_power <- mean(p_values <= alpha)
    mean_time <- mean(computation_times)
    median_time <- median(computation_times)
    
    power_results[i, ] <- c(n, Private_power)
    time_results[i, ] <- c(n, mean_time, median_time)
  }
  
  # Combine results
  result_data <- data.frame(
    n = n_vec,
    Power_or_level = power_results[, "Private_power"],
    Mean_Time = time_results[, "Mean_time"],
    Median_Time = time_results[, "Median_time"]
  )
  
  return(result_data)
}


#Example


# level
joint_probs0 <- c(0.25, 0.25, 0.25, 0.25)  # Null distribution
n_vec <- c(20,30,50, 100, 200, 500, seq(1000, 20000, by = 1000))
epsilon <- 0.1  # Privacy budget
H <- 10^4       # Bootstrap samples
sim <- 10^4     # Monte Carlo simulations
alpha <- 0.05   # Significance level

level_results <- fima_simulation(joint_probs0, n_vec, epsilon, alpha, H, sim)



#power
delta = 0.01  
joint_probs0 <-c(0.25, 0.25,0.25, 0.25) + delta*c(1,-1,-1,1)  # Alternative distribution
n_vec <- c(20,30,50, 100, 200, 500, seq(1000, 20000, by = 1000))
epsilon <- 0.1
H <- 10^4
sim <- 10^4
alpha <- 0.05

power_result <- fima_simulation(joint_probs0, n_vec, epsilon, alpha, H, sim)
