# ==============================================================================
# Power simulation for two-sample proportion tests
# ------------------------------------------------------------------------------


library(dplyr)
library(purrr)
library(progress)
library(tidyverse)
library(poibin)
library(ggplot2)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------
######################################################
#### Helper Functions
######################################################
exact_binom_test <- function(n, theta1, theta2) {
  # Generate samples
  x1 <- rbinom(1, size = n, prob = theta1)
  x2 <- rbinom(1, size = n, prob = theta2)
  
  # Create 2x2 contingency table
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

# Test of test algorithm
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


optimal.m.alpha0 <- function(n, effect, alpha, pub_power, epsilon, theta1, m_grid = NA, ...){
  
  # If no grid for m provided, assign the default grid discussed in the paper
  if(is.na(m_grid[1])){
    m_grid <- c(1:sqrt(n), floor(n/rev(1:(sqrt(n)+1))))
    m_grid <- m_grid[!duplicated(m_grid)]
  }
  
  
  # An function to efficiently compute the power of ToT with balanced sub-samples
  efficient_power <- function(alpha0, m, effect, pub_power, epsilon, n,
                              alpha = 0.05, sub_samp_sizes = NA, ...){
    pub_pow1 <- try(pub_power(n=n,  effect = effect,  alpha = alpha0,p1=theta1))
    if(class(pub_pow1) == "character" | pub_pow1 == 0){ pub_pow1 <- alpha0}
    pub_pow2 <- try(pub_power(n=n, effect = effect,  alpha = alpha0, p1=theta1))
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


pub_power_two_sample_test <- function(n, effect, alpha, p1) {
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


get_data_non_priv = function(theta, u){
  data_set<- (u < theta) *1
  return(data_set)
}


get_pi <- function(theta, eps, n) {
  u <- runif(n)
  w <- runif(1, -0.5, 0.5)
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}


#' FIMA
#'
#' Generates fiducial samples for differentially private inference.
#' 
#' @param pi0 Private proportion estimate
#' @param eps Privacy parameter 
#' @param n Sample size
#' @param H Number of fiducial samples (default = 1)
#' @param delta Sensitivity parameter (default = 1)
#' @param seed Random seed (default = 123)
#' 
#' @return Vector of H fiducial samples
fima <- function(pi0, eps, n, H = 1, delta = 1, seed = 123) {
  set.seed(seed)
  
  # Generate noise variables
  Wj <- runif(H, -0.5, 0.5)
  
  # Create perturbed proportions
  pi_star <- pi0 + delta/(eps*n)*sign(Wj)*log(1 - 2*abs(Wj))
  
  # Initialize fiducial samples
  theta_ifb <- rep(NA, H)
  
  # Generate beta samples for valid proportions
  index <- pi_star < 1 & pi_star > 0
  theta_ifb[index] <- rbeta(
    sum(index), 
    n*pi_star[index] + 0.5,  # Haldane-Anscombe correction
    n*(1 - pi_star[index]) + 0.5
  )
  
  # Handle boundary cases
  theta_ifb[pi_star >= 1] <- 1 - .Machine$double.eps
  theta_ifb[pi_star <= 0] <- .Machine$double.eps
  
  return(theta_ifb)
}


#' Two-Sample FIMA Test
#'
#' Computes the difference between two groups using FIMA.
#' 
#' @param pi1 Private proportion estimate for group 1
#' @param pi2 Private proportion estimate for group 2
#' @param H Number of fiducial samples
#' @param eps Privacy parameter 
#' @param n1 Sample size for group 1
#' @param n2 Sample size for group 2
#' @param seed Random seed
#' 
#' @return Vector of H differences between group distributions
fima2sample <- function(pi1, pi2, H, eps, n1, n2, seed) {
  set.seed(seed)
  
  # Generate distributions for both groups
  distri_grp1 <- fima(pi0 = pi1, H=H, eps = eps/2, n = n1, seed=sample(1:10000,1))
  distri_grp2 <- fima(pi0 = pi2, H=H, eps = eps/2, n = n2, seed=sample(10001:20000,1))
  
  # Compute differences
  res <- distri_grp1 - distri_grp2
  return(res)
}


# 1. Define a unified simulation function
run_combined_simulation <- function(theta1, theta_values, n, B, eps = 1, H = 1000, alpha = 0.05, seed) {
  if (!is.null(seed)) set.seed(seed)
  
  # Initialize results storage
  results <- expand.grid(
    theta = theta_values,
    rep = 1:B,
    method = c("TestOfTest", "NonPrivate", "FIMA", "Exact"),
    rejected = NA,
    time_ms = NA,
    stringsAsFactors = FALSE
  )
  
  pb <- progress_bar$new(total = length(theta_values) * B, format = "[:bar] :percent")
  
  for (j in seq_along(theta_values)) {
    theta <- theta_values[j]
    effect = theta -theta1
    pars <- optimal.m.alpha0(n, effect = effect, alpha, pub_power_two_sample_test, eps, theta, m_grid = NA)
    
    
    for (i in 1:B) {
      set.seed(i)
      
      # --- Test of Test Method ---
      #pars <- optimal.m.alpha0(n, effect = 0, alpha, pub_power_two_sample_test, eps, theta, m_grid = NA)
      data_set1 <- get_data_non_priv(theta1, runif(n))
      data_set <- get_data_non_priv(theta, runif(n))
      subset1 <- data.frame(data_set1)
      subset2 <- data.frame(data_set)
      
      start_time <- Sys.time()
      
      
      test_result <- test.of.tests(
        subset1, subset2, 
        public_two_sampl_prop_test,  
        eps, pars$m, pars$alpha0, 
        sub_samp_sizes = NA
      )
      tot_rejected <- test_result$p.value < alpha
      tot_time <- as.numeric(Sys.time() - start_time) * 1000
      
      # --- Non-Private Method ---
      start_time <- Sys.time()
      np_test <- prop_test(n, theta1 = theta1, theta2 = theta)
      np_rejected <- np_test$p_value < alpha 
      np_time <- as.numeric(Sys.time() - start_time) * 1000
      
      # --- FIMA Method ---
      start_time <- Sys.time()
      pi1 <- get_pi(theta = theta1, eps = eps/2, n = n)
      pi2 <- get_pi(theta = theta, eps = eps/2, n = n)
      fima_samples <- fima2sample(pi1 = pi1, pi2 = pi2, H = H, eps = eps, n1 = n, n2 = n, seed = 10 + i)
      fima_rejected <- mean(fima_samples >= 0) < alpha
      fima_time <- as.numeric(Sys.time() - start_time) * 1000
      
      # --- Exact Method ---
      start_time <- Sys.time()
      exact_test <- exact_binom_test(n, theta1 = theta1, theta2 = theta)
      exact_rejected <- exact_test$p_value < alpha
      exact_time <- as.numeric(Sys.time() - start_time) * 1000
      
      # Store results
      idx <- which(results$theta == theta & results$rep == i)
      results[idx, ] <- data.frame(
        theta = theta,
        rep = i,
        method = c("TestOfTest", "NonPrivate", "FIMA", "Exact"),
        rejected = c(tot_rejected, np_rejected, fima_rejected, exact_rejected),
        time_ms = c(tot_time, np_time, fima_time, exact_time)
      )
      
      pb$tick()
    }
  }
  
  return(results)
}

# Analysis Functions
analyze_results <- function(sim_results) {
  # Theta-specific analysis
  theta_stats <- sim_results %>%
    group_by(theta, method) %>%
    summarise(
      empirical_power = mean(rejected),
      power_se = sd(rejected)/sqrt(n()),
      mean_time_ms = mean(time_ms),
      time_se = sd(time_ms)/sqrt(n()),
      .groups = "drop"
    )  
  return(theta_stats)
}


# Run Simulation
theta_values <- seq(0.2, 0.9, by = 0.1)
theta1 <- 0.2
n <- 30
B <- 10^4  
eps <- 1
H <- 10^4

alpha <- 0.05
m_grid <- NA 
pub_power<-pub_power_two_sample_test


combined_results <- run_combined_simulation(theta1,
  theta_values, n, B, eps, H, seed=123
)


# Create filtered datasets for each method
test_of_test_results <- combined_results %>% filter(method == "TestOfTest")
nonprivate_results <- combined_results %>% filter(method == "NonPrivate")
fima_results <- combined_results %>% filter(method == "FIMA")
exact_results <- combined_results %>% filter(method == "Exact")

# Create analysis function for individual methods
analyze_method <- function(method_data, method_name) {
  # Calculate statistics
  stats <- method_data %>%
    group_by(theta) %>%
    summarise(
      power = mean(rejected),  
      error_se = sd(rejected)/sqrt(n()),
      mean_time_ms = mean(time_ms),
      time_se = sd(time_ms)/sqrt(n()),
      n = n(),
      .groups = "drop"
    )
  return(stats)
  
}

# Analyze each method
tot_analysis <- analyze_method(test_of_test_results, "Test of Test")
np_analysis <- analyze_method(nonprivate_results, "NonPrivate")
fima_analysis <- analyze_method(fima_results, "FIMA")
exact_analysis <- analyze_method(exact_results, "Exact")

tot_analysis
np_analysis
fima_analysis
exact_analysis
