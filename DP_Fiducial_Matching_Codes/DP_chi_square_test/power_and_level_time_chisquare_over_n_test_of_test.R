library(tidyverse)
library(purrr)
library(poibin)
library(progress)



######################################################
#### Helper Functions
######################################################

# Sample from Tulap(m,b)
# Algorithm 2 in Awan & Slavkovic (2018)
rtulap <- function(n, m, b) {
  vect <- rep(NA, n)
  G1 <- rgeom(n, 1 - b)
  G2 <- rgeom(n, 1 - b)
  U <- runif(n, -.5, .5)
  vect <- G1 - G2 + U + m
  return(vect)
}



ptulap <- function(m, b, x){
  if (is.na(x) || is.na(m)) return(NA)
  if (x <= round(m)) {
    (b ^ (- round(x-m))/ (1 + b)) * (b + (x - m - round(x - m) + 1/2) * (1 - b))
  } else {
    1 - ((b ^ round(x - m)) / (1 + b)) * (b + (round(x - m) - (x - m) + 1/2) * (1 - b)) 
  }
}

# Compute p-value given test statistic, Z
# Algorithm 1 in Awan & Slavkovic (2018)
dp.binom.p.val <- function(n, alpha0, epsilon, Z){
  F_underbar <- 0:n %>% map_dbl(~ ptulap(.x, exp(-epsilon), Z))
  B_underbar <- 0:n %>% map_dbl(~ choose(n, .x) * alpha0^.x * (1 - alpha0)^(n - .x))
  return(sum(F_underbar * B_underbar))
}

# Compute probability of statistic as or more extreme than Z under H_A
# where H_A is given by a vector of thetas
dp.binom.alt.prob <- function(n, thetas, epsilon, Z){
  F_underbar <- 0:n %>% map_dbl(~ ptulap(.x, exp(-epsilon), Z))
  B_underbar <- 0:n %>% map_dbl(.f = dpoibin, pp = thetas)
  return(sum(F_underbar * B_underbar))
}




# Public Chi-squared test for independence
public_test <- function(table) {
  
  # Perform chi-square test of independence
  test_result <- chisq.test(table, correct = FALSE)
  p_value <- test_result$p.value  # Store the p-valu
  
  
  return(p_value)
}




pub_power_chisq <- function(n, effect_size, alpha,r=2, c=2) {
  # Degrees of freedom
  df <- (r - 1) * (c - 1)
  
  # Critical value from chi-squared distribution
  critical_value <- qchisq(1 - alpha, df)
  
  # Compute the non-centrality parameter (lambda)
  lambda <- n * effect_size
  
  # Compute the power: 1 - Type II error probability
  power <- 1 - pchisq(critical_value, df, ncp = lambda)
  
  return(power)
}



######## Data Generation  ############
data_generation<-function(n, joint_probs){
  outcomes <- sample(1:4, n, replace = TRUE, prob = joint_probs)
  return(outcomes)
}


optimal_m_alpha0 <- function(effect, pub_power, epsilon, n, m_grid = NA, 
                             alpha = 0.05, ...){
  # If no grid for m provided, assign the default grid discussed in the paper
  if(is.na(m_grid[1])){
    m_grid <- c(1:sqrt(n), floor(n/rev(1:(sqrt(n)+1))))
    m_grid <- m_grid[!duplicated(m_grid)]
  }
  
  # An function to efficiently compute the power of ToT with balanced sub-samples
  efficient_power <- function(alpha0, m, effect, pub_power, epsilon, n,  
                              alpha = 0.05, sub_samp_sizes = NA, ...){
    pub_pow1 <- try(pub_power(effect = effect, n = ceiling(n/m), alpha = alpha0, ...))
    if(class(pub_pow1) == "character" | pub_pow1 == 0){ pub_pow1 <- alpha0}
    pub_pow2 <- try(pub_power(effect = effect, n = floor(n/m), alpha = alpha0, ...))
    if(class(pub_pow2) == "character" | pub_pow2 == 0){ pub_pow2 <- alpha0}
    thetas <- c(rep(pub_pow1, floor(n - floor(n/m)*m)),
                rep(pub_pow2, ceiling(m - n + floor(n/m)*m)))
    tol <- 2*log(100)/epsilon
    critical_value <- optimize(function(x){abs(dp.binom.p.val(m, alpha0, epsilon, x) - (1 - alpha))},
                               interval = c(-tol, m+tol), tol = 5e-4)$minimum
    return(1 - dp.binom.alt.prob(m, thetas, epsilon, critical_value))
  }
  
  # Function to find the alpha0 that maximizes power for a given m
  optimal_alpha0 <- function(m, n, effect, epsilon, alpha, pub_power, ...){
    opt <- optimize(f = efficient_power, interval = c(0,1), maximum = T, tol = 5e-3,
                    m = m, effect = effect, pub_power = pub_power, epsilon = epsilon,
                    n = n, alpha = alpha, ...)
    return(c(opt$objective, opt$maximum, m))
  }
  
  # For each m in grid, find the alpha0 that maximizes power
  x <- lapply(X = m_grid, FUN = optimal_alpha0, n = n, effect = effect, 
              epsilon = epsilon, alpha = alpha, pub_power = pub_power, ...)
  
  # Find the m,alpha0 combination that yields the maximum overall power and return summary
  max_x <- x[[which.max(as.numeric(sapply(x,"[[",1)))]]
  return(list("power" = max_x[1], "alpha0" = max_x[2], "m" = as.integer(max_x[3])))
}




practical_m_alpha0 <- function(rho, pub_power, epsilon, n, m_grid = NA, 
                               effect_grid = NA, alpha = 0.05, ...){
  # If no grid for m provided, assign the default grid discussed in the paper
  if(is.na(m_grid[1])){
    m_grid <- c(1:sqrt(n), floor(n/rev(1:(sqrt(n)+1))))
    m_grid <- m_grid[!duplicated(m_grid)]
  }
  # If no grid for effect provided, assign a default grid
  if(is.na(effect_grid[1])){
    effect_grid <- 2^c(-7:-1, seq(0,4,0.5))
  }
  
  # Perform a binary search for the smallest effect that gives power greater than rho
  i <- (length(effect_grid)+1) %/% 2; power <- 0; effect <- Inf
  while(length(effect_grid) > 1){
    opt <- optimal_m_alpha0(effect = effect_grid[i], pub_power = pub_power, 
                            epsilon = epsilon, n = n, m_grid = m_grid, 
                            alpha = alpha, ...)
    if(opt$power >= rho){
      power <- opt$power; alpha0 <- opt$alpha0; m <- opt$m; effect <- effect_grid[i]
      effect_grid <- effect_grid[1:i]
      i <- (i+1) %/% 2
      
    }
    else{
      effect_grid <- effect_grid[(i+1):length(effect_grid)]
      i <- (length(effect_grid) + 1) %/% 2
    }
  }
  if(effect != effect_grid){
    opt <- optimal_m_alpha0(effect = effect_grid, pub_power = pub_power, 
                            epsilon = epsilon, n = n, m_grid = m_grid, 
                            alpha = alpha, ...)
    power <- opt$power; alpha0 <- opt$alpha0; m <- opt$m
  }
  return(list("power" = power, "alpha0" = alpha0, "m" = m, "effect" = effect_grid))
}


############ Test of Test Algorithm #################

test_of_tests <- function(outcomes, n, tau, epsilon, alpha, m, alpha0) {
  
  # Partition outcomes
  subset_size <- floor(n / m)
  subset_indices <- split(outcomes[1:(m * subset_size)], rep(1:m, each = subset_size))
  if (n %% m != 0) {
    subset_indices[[m]] <- c(subset_indices[[m]], outcomes[(m * subset_size + 1):n])
  }
  
  # Generate contingency tables
  tables <- lapply(subset_indices, function(subset_outcomes) {
    table <- matrix(0, nrow = 2, ncol = 2)
    for (outcome in subset_outcomes) {
      row_index <- (outcome - 1) %/% 2 + 1
      col_index <- (outcome - 1) %% 2 + 1
      table[row_index, col_index] <- table[row_index, col_index] + 1
    }
    return(table)
  })
  
  # Compute test statistics
  p_values <- suppressWarnings(
    
    sapply(tables, function(table) {
      if(tau(table)=="character"){
        runif(1)
        
      }
      else{
        tau(table)
      }
    })
  )
  
  # Private binomial test
  Z <- rtulap(n = 1, m = sum(p_values <= alpha0), b = exp(-epsilon))
  p_val <- 1 - dp.binom.p.val(n = m, alpha0 = alpha0, epsilon = epsilon, Z = Z)
  
  return(list("Z" = Z, "p_value" = p_val))
}





simulation_test_of_test_chi_squared <- function(joint_probs_H_0, joint_probs_H_a, rho, n_sim, n_values, 
                                                pub_test, pub_power, epsilon, alpha, m_grid, effect_grid){
  
  
  
  # Initialize matrices to store results
  result_matrix_H_0 <- matrix(NA, nrow = n_sim, ncol = length(n_values))
  result_matrix_H_a <- matrix(NA, nrow = n_sim, ncol = length(n_values))
  time_matrix_H_0 <- matrix(NA, nrow = n_sim, ncol = length(n_values))
  time_matrix_H_a <- matrix(NA, nrow = n_sim, ncol = length(n_values))
  
  # Create progress bar
  total_iterations <- length(n_values) * n_sim
  pb <- progress_bar$new(
    format = "  n=:n_val [:bar] :percent | :current/:total sims | Elapsed: :elapsed | ETA: :eta",
    total = total_iterations,
    clear = FALSE,
    width = 80,
    show_after = 0
  )
  
  for (i in seq_along(n_values)) {
    n <- n_values[i]
    pars <- tryCatch({
      practical_m_alpha0(rho, pub_power, epsilon, n, m_grid, effect_grid, alpha)
    }, error = function(e) {
      message("\nError in practical_m_alpha0 for n=", n, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(pars)) next
    
    m <- pars$m
    alpha0 <- pars$alpha0
    
    for (j in 1:n_sim) {
      set.seed(j * 10 + j)
      pb$tick(tokens = list(n_val = n))
      
      tryCatch({
        # Type I error calculation
        start_time_H_0 <- Sys.time()
        outcomes_H_0 <- data_generation(n, joint_probs_H_0)
        result_matrix_H_0[j,i] <- test_of_tests(outcomes_H_0, n, tau=pub_test, epsilon, alpha, m, alpha0)$p_value
        time_matrix_H_0[j,i] <- as.numeric(Sys.time() - start_time_H_0) * 1000
        
        # Power calculation
        start_time_H_a <- Sys.time()
        outcomes_H_a <- data_generation(n, joint_probs_H_a)
        result_matrix_H_a[j,i] <- test_of_tests(outcomes_H_a, n, tau=pub_test, epsilon, alpha, m, alpha0)$p_value
        time_matrix_H_a[j,i] <- as.numeric(Sys.time() - start_time_H_a) * 1000
      }, error = function(e) {
        message("\nError in simulation ", j, " for n=", n, ": ", e$message)
      })
    }
  }
  
  # Calculate results
  Type_I_error <- sapply(1:length(n_values), function(i) {
    vals <- result_matrix_H_0[,i]
    mean(vals < alpha, na.rm = TRUE)
  })
  
  power <- sapply(1:length(n_values), function(i) {
    vals <- result_matrix_H_a[,i]
    mean(vals < alpha, na.rm = TRUE)
  })
  
  # Create results dataframe
  results_df <- data.frame(
    n = n_values,
    power = ifelse(is.nan(power), NA, power),
    Type_I_error = ifelse(is.nan(Type_I_error), NA, Type_I_error),
    mean_time_power_ms = apply(time_matrix_H_a, 2, mean, na.rm = TRUE),
    median_time_power_ms = apply(time_matrix_H_a, 2, median, na.rm = TRUE),
    mean_time_typeI_ms = apply(time_matrix_H_0, 2, mean, na.rm = TRUE),
    median_time_typeI_ms = apply(time_matrix_H_0, 2, median, na.rm = TRUE)
  )
  
  # Replace NA/NaN with 0 for timing columns
  results_df[is.na(results_df)] <- 0
  
  #############################################
  # PLOTTING SECTION 
  #############################################
  
  # Plot 1: Power with computation times
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 4) + 0.3)
  
  # Power plot
  ylim_power <- range(c(0, 1, power), finite = TRUE)
  ylim_time_power <- range(c(0, results_df$mean_time_power_ms, results_df$median_time_power_ms), finite = TRUE)
  
  plot(n_values, power, type = "b", pch = 16, col = "blue", 
       xlab = "Sample Size (n)", ylab = "Power", 
       main = "Power and Computation Time",
       ylim = ylim_power)
  par(new = TRUE)
  plot(n_values, results_df$mean_time_power_ms, type = "b", pch = 17, col = "red", 
       axes = FALSE, xlab = "", ylab = "", 
       ylim = ylim_time_power)
  par(new = TRUE)
  plot(n_values, results_df$median_time_power_ms, type = "b", pch = 18, col = "green", 
       axes = FALSE, xlab = "", ylab = "", 
       ylim = ylim_time_power)
  axis(side = 4)
  mtext("Time (ms)", side = 4, line = 3)
  legend("topleft", legend = c("Power", "Mean Time", "Median Time"),
         col = c("blue", "red", "green"), pch = c(16, 17, 18), lty = 1)
  
  # Plot 2: Type I Error with computation times
  ylim_typeI <- range(c(0, max(Type_I_error, na.rm = TRUE)), finite = TRUE)
  ylim_time_typeI <- range(c(0, results_df$mean_time_typeI_ms, results_df$median_time_typeI_ms), finite = TRUE)
  
  plot(n_values, Type_I_error, type = "b", pch = 16, col = "blue", 
       xlab = "Sample Size (n)", ylab = "Type I Error", 
       main = "Type I Error and Computation Time",
       ylim = ylim_typeI)
  par(new = TRUE)
  plot(n_values, results_df$mean_time_typeI_ms, type = "b", pch = 17, col = "red", 
       axes = FALSE, xlab = "", ylab = "", 
       ylim = ylim_time_typeI)
  par(new = TRUE)
  plot(n_values, results_df$median_time_typeI_ms, type = "b", pch = 18, col = "green", 
       axes = FALSE, xlab = "", ylab = "", 
       ylim = ylim_time_typeI)
  axis(side = 4)
  mtext("Time (ms)", side = 4, line = 3)
  legend("topleft", legend = c("Type I Error", "Mean Time", "Median Time"),
         col = c("blue", "red", "green"), pch = c(16, 17, 18), lty = 1)
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
  
  return(list(
    results = results_df,
    Type_I_error = Type_I_error,
    power = power,
    time_matrix_H_0 = time_matrix_H_0,
    time_matrix_H_a = time_matrix_H_a
  ))
}



# Parameters
delta = 0.01
n_values <- c(20,30,50, 100, 200, 500, seq(1000, 20000, by = 1000))
joint_probs_H_0 = c(0.25, 0.25, 0.25, 0.25)
joint_probs_H_a = c(0.25, 0.25, 0.25, 0.25) + delta * c(1, -1, -1, 1) 
rho = 0.90
n_sim = 10^3  
epsilon = 0.5
alpha = 0.05
m_grid = 1:50
effect_grid <- seq(0.01, 3, 0.01)
pub_test = public_test
pub_power = pub_power_chisq



# Run simulation with error handling
tryCatch({
  power_and_level_test_of_test <- simulation_test_of_test_chi_squared(
    joint_probs_H_0, joint_probs_H_a, rho, n_sim, n_values, 
    pub_test, pub_power, epsilon, alpha, m_grid, effect_grid
  )
  print(power_and_level_test_of_test$results)
}, error = function(e) {
  message("Error occurred: ", e$message)
  traceback()
})


