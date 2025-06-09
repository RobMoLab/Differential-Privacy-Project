library(tidyverse)
library(purrr)
library(poibin)


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

# CDF of Tulap(m,b)
# Definition 4.1 in Awan & Slavkovic (2018)
ptulap <- function(m, b, x){
  if (x <= round(m)) {
    (b ^ (- round(x-m))/ (1 + b)) * (b + (x - m - round(x - m) + 1/2) * (1 - b))
  } else {
    1 - ((b ^ round(x - m)) / (1 + b)) * (b + (round(x - m) - (x - m) + 1/2) * (1 - b))  }
}

# Compute p-value given test statistic, Z
# Algorithm 1 in Awan & Slavkovic (2018)
dp.binom.p.val <- function(n, alpha0, epsilon, Z){
  F_underbar <- 0:n %>%
    map_dbl(~ ptulap(.x, exp(-epsilon), Z))
  B_underbar <- 0:n %>%
    map_dbl(~ choose(n, .x) * alpha0^.x * (1 - alpha0)^(n - .x))
  return(sum(F_underbar * B_underbar))
}

# Compute probability of statistic as or more extreme than Z under H_A
# where H_A is given by a vector of thetas
dp.binom.alt.prob <- function(n, thetas, epsilon, Z){
  F_underbar <- 0:n %>%
    map_dbl(~ ptulap(.x, exp(-epsilon), Z))
  B_underbar <- 0:n %>%
    map_dbl(.f = dpoibin, pp = thetas)
  return(sum(F_underbar * B_underbar))
}




test.of.tests <- function(subset1, subset2, pub_test, epsilon, m, alpha0, sub_samp_sizes = NA) {
  # Ensure subsets are not NULL or empty
  if (is.null(subset1) || is.null(subset2)) stop("subset1 and subset2 cannot be NULL.")
  if (length(subset1) == 0 || length(subset2) == 0) stop("subset1 and subset2 cannot be empty.")
  
  # Handle vectors and matrices appropriately
  n <- if (is.matrix(subset1) || is.data.frame(subset1)) nrow(subset1) else length(subset1)
  
  # Ensure valid n and m
  if (is.na(n) || is.na(m) || m <= 0 || n <= 0) stop("n and m must be positive integers.")
  
  if (m > n) {
    warning("Reducing m to match n because m > n.")
    m <- n
  }
  # m=6
  # n=100
  if (!is.na(sub_samp_sizes[1])) {
    sub_samples <- sample(rep(1:m, sub_samp_sizes))
  } else {
    sub_samples <- sample(rep(1:m, length.out = n))
  }
  
  # Initialize p-value vector
  sub_tests <- rep(NA, m)
  
  for (i in 1:m) {
    # Compute p-value in each subsample
    sub_tests[i] <- tryCatch(
      pub_test(subset1[sub_samples == i, ], subset2[sub_samples == i, ]),
      error = function(e) runif(1)
    )
    
    # Ensure result is numeric; if not, sample from Unif(0, 1)
    if (!is.numeric(sub_tests[i]) || is.na(sub_tests[i])) {
      sub_tests[i] <- runif(1)
    }
  }
  
  # Run the Awan & Slavkovic (2018) private binomial test
  Z <- rtulap(n = 1, m = sum(sub_tests <= alpha0), b = exp(-epsilon))
  p_val <- 1 - dp.binom.p.val(n = m, alpha0 = alpha0, epsilon = epsilon, Z = Z)
  
  return(list("Z" = Z, "p.value" = p_val))
}








optimal.m.alpha0 <- function(n, effect, alpha, pub_power, epsilon, m_grid = NA, ...){
  
  # If no grid for m provided, assign the default grid discussed in the paper
  if(is.na(m_grid[1])){
    m_grid <- c(1:sqrt(n), floor(n/rev(1:(sqrt(n)+1))))
    m_grid <- m_grid[!duplicated(m_grid)]
  }
  
  
  # An function to efficiently compute the power of ToT with balanced sub-samples
  efficient_power <- function(alpha0, m, effect, pub_power, epsilon, n,
                              alpha = 0.05, sub_samp_sizes = NA, ...){
    pub_pow1 <- try(pub_power(n=n,  effect = effect,  alpha = alpha0))
    if(class(pub_pow1) == "character" | pub_pow1 == 0){ pub_pow1 <- alpha0}
    pub_pow2 <- try(pub_power(n=n, effect = effect,  alpha = alpha0))
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
              epsilon = epsilon, alpha = alpha, pub_power = pub_power)
  
  # Find the m,alpha0 combination that yields the maximum overall power and return summary
  max_x <- x[[which.max(as.numeric(sapply(x,"[[",1)))]]
  return(list("power" = max_x[1], "alpha0" = max_x[2], "m" = as.integer(max_x[3])))
}




practical.m.alpha0 <- function(rho, n, alpha, pub_power, epsilon, m_grid = NA,
                               effect_grid = NA, ...){
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
    opt <- optimal.m.alpha0( n=n, effect = effect_grid[i], alpha = alpha, pub_power = pub_power,
                             epsilon = epsilon, m_grid = m_grid)
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
    opt <- optimal.m.alpha0(n=n, effect = effect_grid, alpha = alpha, pub_power = pub_power,
                            epsilon = epsilon, m_grid = m_grid)
    power <- opt$power; alpha0 <- opt$alpha0; m <- opt$m
  }
  return(list("power" = power, "alpha0" = alpha0, "m" = m, "effect" = effect_grid))
}







public_two_sampl_prop_test <- function(subset1, subset2) {
  n1 <- length(subset1)
  n2 <- length(subset2)
  x1<-sum(subset1)
  x2<-sum(subset2)
  
  # Compute sample proportions
  p1 <- x1 / n1
  p2 <- x2 / n2
  
  # Null hypothesis: H0: theta1 - theta2 = 0
  # Alternative hypothesis: H1: theta1 - theta2 < 0
  
  # Calculate pooled proportion under H0
  p_pooled <- (x1 + x2) / (n1 + n2)
  
  # Compute the test statistic (z-score)
  z <- (p1 - p2) / sqrt(p_pooled * (1 - p_pooled) * (1/n1 + 1/n2))
  
  # Compute the p-value for the one-tailed test
  p_value <- pnorm(z)
  
  return(p_value)
}




pub_power_two_sample_test <- function(n, effect, alpha, p1 = 0.8) {
  p2 <- effect + p1
  
  # Clip probabilities to [0,1] to avoid NaNs
  p1 <- max(0, min(1, p1))
  p2 <- max(0, min(1, p2))
  
  # Compute power (return alpha0 if power.prop.test fails)
  tryCatch({
    power <- power.prop.test(
      n = n,
      p1 = p1,
      p2 = p2,
      sig.level = alpha,
      alternative = "one.sided"
    )$power
    if (is.na(power)) return(alpha)  # Fallback to alpha if power=NA
    return(power)
  }, error = function(e) {
    return(alpha)  # Return alpha if error (e.g., p1=p2)
  })
}












simulation_test_of_test <- function(n_sim, theta1, theta2, n_values, alpha, epsilon, pub_test, pub_power, m_grid, seed) {
  # Initialize matrices to store results
  p_value_H_a <- matrix(NA, nrow = n_sim, ncol = length(n_values))
  time_matrix <- matrix(NA, nrow = n_sim, ncol = length(n_values))  # To store computation times
  
  set.seed(seed)
  effect <- abs(theta2 - theta1)
  
  for (j in seq_along(n_values)) {
    n <- n_values[j]
    
    # Get parameters once per sample size (outside the simulation loop)
    pars <- optimal.m.alpha0(n, effect, alpha = alpha, pub_power = pub_power, 
                             epsilon = epsilon, m_grid = m_grid)
    m <- pars$m
    alpha0 <- pars$alpha0
    
    for (i in 1:n_sim) {
      # Start timing
      start_time <- Sys.time()
      
      # Generate data under the alternative hypothesis
      subset1_H_a <- data.frame(rbinom(n, 1, theta1))
      subset2_H_a <- data.frame(rbinom(n, 1, theta2))
      
      # Perform the test
      test_result <- test.of.tests(subset1_H_a, subset2_H_a, pub_test, epsilon, m, alpha0, sub_samp_sizes = NA)
      p_value_H_a[i, j] <- test_result$p.value
      
      # Record time in milliseconds
      time_matrix[i, j] <- as.numeric(Sys.time() - start_time) * 1000
    }
  }
  
  # Calculate results
  power <- apply(p_value_H_a < alpha, 2, mean)
  mean_times <- apply(time_matrix, 2, mean)
  median_times <- apply(time_matrix, 2, median)
  
  # Create results dataframe
  results_df <- data.frame(
    n = n_values,
    power = power,
    mean_time_ms = mean_times,
    median_time_ms = median_times
  )
  
  # Create the combined plot
  par(mar = c(5, 4, 4, 4) + 0.3)  
  
  # Plot power (left y-axis)
  plot(n_values, power, type = "b", pch = 16, col = "blue", 
       xlab = "Sample Size (n)", ylab = "Power", 
       main = "Power and Computation Time vs Sample Size",
       ylim = c(0, 1))
  
  # Add mean times (right y-axis)
  par(new = TRUE)
  plot(n_values, mean_times, type = "b", pch = 17, col = "red", 
       axes = FALSE, xlab = "", ylab = "", ylim = c(0, max(mean_times, median_times) * 1.1))
  
  # Add median times (right y-axis)
  par(new = TRUE)
  plot(n_values, median_times, type = "b", pch = 18, col = "green", 
       axes = FALSE, xlab = "", ylab = "", ylim = c(0, max(mean_times, median_times) * 1.1))
  
  # Add right y-axis
  axis(side = 4, at = pretty(range(c(mean_times, median_times))))
  mtext("Computation Time (ms)", side = 4, line = 3)
  
  # Add legend
  legend("topleft", legend = c("Power", "Mean Time", "Median Time"), 
         col = c("blue", "red", "green"), pch = c(16, 17, 18), lty = 1)
  
  return(list(
    results = results_df,
    power = power,
    mean_times = mean_times,
    median_times = median_times
  ))
}

# Parameters
theta1 = 0.8
theta2 = 0.9
epsilon <- 1
alpha <- 0.05
m_grid <- NA
seed <- 134
n_values <- c(16, 30, 50, 100, 150, 200, 350, 400, 500)
pub_test = public_two_sampl_prop_test
pub_power = pub_power_two_sample_test
n_sim <- 10^4

# Run simulation
sim_results <- simulation_test_of_test(n_sim, theta1, theta2, n_values, alpha, epsilon, 
                                       pub_test, pub_power, m_grid, seed)


#power_and_time_two_samples_over_n_test_of_test<-sim_results
# Save results
#save(power_and_time_two_samples_over_n_test_of_test, file = "power_and_time_two_samples_over_n_test_of_test.Rda")


# Print the results dataframe
print(sim_results$results)

# The plot is automatically generated by the function
