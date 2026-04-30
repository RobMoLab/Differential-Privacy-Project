# ==============================================================================
# Power and level simulation for DP chi-square tests
# ------------------------------------------------------------------------------


library(glmnet)
library(nloptr)
library(alabama)
library(parallel)
library(pbmcapply)  
library(tidyverse)
library(purrr)
library(poibin)
library(progress)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------
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


################ Gaboardi Method     ##########################


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


##### Fima method   #########
#' Compute Chi-Squared Statistic
#'
#' Calculates the chi-squared statistic for a contingency table.
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


#' Generate Laplace Noise
#'
#' @param sensitivity The sensitivity of the query
#' @param epsilon Privacy parameter
#' @return Laplace noise value
laplace_noise <- function(sensitivity, epsilon){
  W <- runif(1, -0.5, 0.5)
  noise <- -sensitivity/epsilon*sign(W)*log(1 - 2*abs(W))
  return(noise)
}

#' Apply Laplace Mechanism to Count Data
#'
#' @param counts Vector of counts to privatize
#' @param epsilon Privacy parameter
#' @return Noisy counts with non-negative values
add_laplace_noise <- function(counts, epsilon) {
  sensitivity <- 2
  noisy_counts <- sapply(counts, function(count) count+laplace_noise(sensitivity, epsilon=epsilon))
  noisy_counts <- pmax(noisy_counts, 0)
  return(haldane_anscombe_correction(noisy_counts))
}

#' Haldane-Anscombe Correction
#'
#' Adjusts zero counts to avoid computational issues.
#' @param counts Vector of counts
#' @return Adjusted counts
haldane_anscombe_correction <- function(counts){
  counts[counts==0] <- counts[counts==0]+0.5
  return(counts)
}

#' Private Proportion Estimator
#'
#' @param prop Original proportion
#' @param n Sample size
#' @param epsilon Privacy parameter
#' @return Noisy proportion
priv_prop <- function(prop, n, epsilon){
  W <- runif(1, -0.5, 0.5)
  noisy_prop <- prop - 2/(n*epsilon)*sign(W)*log(1 - 2*abs(W))
  return(noisy_prop)
}


########## Test of test method ###########
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


#### Simulation #######

combined_simulation <- function(joint_probs0, n_vec, epsilon, delta, alpha, H, k, sim, 
                                rho = 0.9, m_grid = 1:50, 
                                effect_grid = seq(0.01, 3, 0.01)) {
  
  # Initialize results matrices
  results <- matrix(NA, nrow = length(n_vec), ncol = 9)
  colnames(results) <- c("Sample_size", "FIMA_power", "ToT_power", "Gabroadi_power", "Nonprivate_power",
                         "FIMA_time", "ToT_time", "Gabroadi_time", "Nonprivate_time")
  
  # For time matrices (detailed timing info)
  time_matrices <- list(
    FIMA = matrix(NA, nrow = sim, ncol = length(n_vec)),
    ToT = matrix(NA, nrow = sim, ncol = length(n_vec)),
    Gaboardi = matrix(NA, nrow = sim, ncol = length(n_vec)),
    Nonprivate = matrix(NA, nrow = sim, ncol = length(n_vec))
  )
  
  # Create progress bar
  pb <- progress_bar$new(
    format = "  n=:n_val [:bar] :percent | :current/:total sims | Elapsed: :elapsed | ETA: :eta",
    total = length(n_vec) * sim,
    clear = FALSE,
    width = 80,
    show_after = 0
  )
  
  for (i in seq_along(n_vec)) {
    n <- n_vec[i]
    
    # Get optimal parameters for Test of Tests method
    pars <- practical_m_alpha0(rho, pub_power_chisq, epsilon, n, m_grid, effect_grid, alpha)
    m <- pars$m
    alpha0 <- pars$alpha0
    
    # Parallel computation
    sim_results <- pbmclapply(1:sim, function(b) {
      set.seed(b * 10 + b)
      
      # Initialize result storage
      res <- list(
        fima = c(NA, NA),
        tot = c(NA, NA),
        gab = c(NA, NA),
        np = c(NA, NA)
      )
      
      tryCatch({
        # Generate data once to be used by all methods
        counts <- rmultinom(1, n, prob = joint_probs0)
        contingency_table <- matrix(counts, nrow = 2)
        
        ### FIMA Method ###
        start_time <- Sys.time()
        noisy_contingency_table <- matrix(add_laplace_noise(counts, epsilon), nrow = 2)
        chi2_0b <- compute_chi2(noisy_contingency_table)
        
        # Bootstrap under null
        denoised_counts <- noisy_contingency_table
        n_star <- sum(denoised_counts)
        marginal_row <- apply(denoised_counts, 1, sum) / n_star
        marginal_col <- apply(denoised_counts, 2, sum) / n_star
        
        # Generate fiducial samples
        marginal_row_tilde <- matrix(rbeta(H*length(marginal_row), 
                                           n_star*marginal_row + 0.5, 
                                           n_star*(1 - marginal_row) + 0.5),
                                     nrow = H, byrow = TRUE)
        marginal_col_tilde <- matrix(rbeta(H*length(marginal_col), 
                                           n_star*marginal_col + 0.5, 
                                           n_star*(1 - marginal_col) + 0.5),
                                     nrow = H, byrow = TRUE)
        
        # Compute bootstrap samples
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
        
        res$fima <- c(mean(chi2_hb >= chi2_0b), as.numeric(Sys.time() - start_time) * 1000)
        
        ### Test of Tests Method ###
        start_time <- Sys.time()
        outcomes <- data_generation(n, joint_probs0)
        tot_result <- test_of_tests(outcomes, n, public_test, epsilon, alpha, m, alpha0)
        res$tot <- c(tot_result$p_value, as.numeric(Sys.time() - start_time) * 1000)
        
        ### Gabroadi Method ###
        start_time <- Sys.time()
        gab_result <- tryCatch({
          mc_indep_test(contingency_table, epsilon, delta, alpha, "Laplace", k)
        }, error = function(e) {
          "Fail to Reject H0: The variables are independent."
        })
        res$gab <- c(ifelse(gab_result == "Reject H0: The variables are dependent.", 1, 0),
                     as.numeric(Sys.time() - start_time) * 1000)
        
        ### Non-private Method ###
        start_time <- Sys.time()
        np_test <- chisq.test(contingency_table, correct = FALSE)
        res$np <- c(np_test$p.value, as.numeric(Sys.time() - start_time) * 1000)
        
      }, error = function(e) {
        message("Error in iteration ", b, ": ", e$message)
      })
      
      pb$tick(tokens = list(n_val = n))
      return(res)
    }, mc.cores = detectCores())
    
    # Process results with error handling
    process_results <- function(method) {
      vals <- sapply(sim_results, function(x) {
        if (is.null(x[[method]])) return(c(NA, NA))
        x[[method]]
      })
      list(
        pvals = vals[1,],
        times = vals[2,]
      )
    }
    
    fima <- process_results("fima")
    tot <- process_results("tot")
    gab <- process_results("gab")
    np <- process_results("np")
    
    # Store results
    results[i, ] <- c(
      n,
      mean(fima$pvals <= alpha, na.rm = TRUE),  # FIMA power
      mean(tot$pvals <= alpha, na.rm = TRUE),   # ToT power
      mean(gab$pvals, na.rm = TRUE),            # Gabroadi power
      mean(np$pvals <= alpha, na.rm = TRUE),    # Non-private power
      mean(fima$times, na.rm = TRUE),           # FIMA mean time
      mean(tot$times, na.rm = TRUE),            # ToT mean time
      mean(gab$times, na.rm = TRUE),            # Gabroadi mean time
      mean(np$times, na.rm = TRUE)              # Non-private mean time
    )
    
    # Store detailed timing
    time_matrices$FIMA[,i] <- fima$times
    time_matrices$ToT[,i] <- tot$times
    time_matrices$Gaboardi[,i] <- gab$times
    time_matrices$Nonprivate[,i] <- np$times
  }
  
  return(list(
    summary_results = as.data.frame(results),
    time_matrices = time_matrices
  ))
}

# Parameters for simulation
delta <- 0.01
joint_probs_H0 <- c(0.25, 0.25, 0.25, 0.25)
joint_probs_Ha <- c(0.25, 0.25, 0.25, 0.25) + delta * c(1, -1, -1, 1)
n_vec <- c(20, 30, 50, 100, 200, 500, seq(1000, 20000, by = 1000))
epsilon <- 0.5
delta_dp <- 1e-5
alpha <- 0.05
H <- 10^4  # Bootstrap samples for FIMA
k <- 10^3  # MC samples for Gabroadi
sim <- 10^3  # number of simulations

# Run null simulation (Type I error)
level_results <- combined_simulation(joint_probs_H0, n_vec, epsilon, delta_dp, alpha, H, k, sim)
level_results$summary_results
# Run alternative simulation (Power)
power_results <- combined_simulation(joint_probs_Ha, n_vec, epsilon, delta_dp, alpha, H, k, sim)

power_results$summary_results
