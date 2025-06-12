# Load necessary package
library(microbenchmark)
library(parallel)

# Function Definitions
inverse_logit <- function(nu) {
  exp(nu) / (1 + exp(nu))
}




get_pi <- function(p_hat, eps, n) {
  w = runif(1, -0.5, 0.5)
  # Ensure w is strictly between -0.5 and 0.5
  if (w == -0.5) w <- -0.5 + .Machine$double.eps
  if (w == 0.5) w <- 0.5 - .Machine$double.eps
  p_hat - 1 / (eps * n) * sign(w) * log(1 - 2 * abs(w))
}


JINI_algorithm_Fiducial <- function(pi0, H, eps, n, err=0.01, seed) {
  set.seed(seed)
  res <- rep(NA, H)
  counter <- 0
  while (counter < H) {
    if(!is.na(eps)){
      Wj <- runif(1, -0.5, 0.5)
      deltaj <- pi0 + 1 / (eps * n) * sign(Wj) * log(1 - 2 * abs(Wj))
    }else{
      deltaj <- pi0
    }
    if (!is.na(deltaj)) {
      if(deltaj>=1) deltaj = 1-err
      if(deltaj<=0) deltaj = err
      counter <- counter + 1
      thetaj <- rbeta(1, n * deltaj + 0.5, n * (1 - deltaj) + 0.5)
      res[counter] <- thetaj
    }
  }
  res
}


# Main function with parallelization
main_function <- function() {
  # Parameters
  eps_values <- c(0.1, 0.3, 1, 3, 10, Inf)
  N <- c(15, 30, 100, 200, 500, 1000, 2000)
  alpha <- 0.05
  H <- 10^4
  beta0 <- 0.5
  beta1 <- 2
  beta <- c(beta0, beta1)
  p1 <- 0.5
  
  # Results Storage
  initialize_matrices <- function(dim1, dim2) {
    matrix(NA, dim1, dim2)
  }
  
  beta0_coverage <- initialize_matrices(length(N), length(eps_values))
  beta1_coverage <- initialize_matrices(length(N), length(eps_values))
  beta0_CI_length <- initialize_matrices(length(N), length(eps_values))
  beta1_CI_length <- initialize_matrices(length(N), length(eps_values))
  
  # Main Loops
  for (k in 1:length(eps_values)) {
    for (m in 1:length(N)) {
      set.seed(m + 112)
      X1 <- rbinom(N[m], 1, p1)
      nu <- beta0 + beta1 * X1
      probs <- inverse_logit(nu)
      seed <- N[m] + 112 
      
      # Simulation with Parallelization
      B <- 10^4
      
      simulation_results <- mclapply(1:B, function(i) {
        set.seed(i)
        u <- runif(N[m])
        Y <- as.numeric(probs > u)
        
        # Private JINI Fiducial
        pi1 <- mean(Y[X1 == 0])
        pi2 <- mean(Y[X1 == 1])
        
        n1 <- sum(X1 == 0)
        n2 <- sum(X1 == 1)
        
        pi1_priv <- get_pi(pi1, eps_values[k], n1)
        pi2_priv <- get_pi(pi2, eps_values[k], n2)
        
        jini_dis_pi1 <- JINI_algorithm_Fiducial(pi1_priv, H, eps_values[k], n1, seed=seed + 2 * i)
        jini_dis_pi2 <- JINI_algorithm_Fiducial(pi2_priv, H, eps_values[k], n2, seed=seed + 2 * i)
        
        beta_0 <- log(jini_dis_pi1 / (1 - jini_dis_pi1))
        beta_1 <- log(jini_dis_pi2 / (1 - jini_dis_pi2)) - beta_0
        
        CI_beta0 <- quantile(beta_0, c(alpha / 2, 1 - alpha / 2))
        CI_beta1 <- quantile(beta_1, c(alpha / 2, 1 - alpha / 2))
        
        len_CI_beta0 <- CI_beta0[2] - CI_beta0[1]
        len_CI_beta1 <- CI_beta1[2] - CI_beta1[1]
        
        cov_beta0 <- (CI_beta0[1] < beta0) & (CI_beta0[2] > beta0)
        cov_beta1 <- (CI_beta1[1] < beta1) & (CI_beta1[2] > beta1)
        
        list(cov_beta0 = cov_beta0, cov_beta1 = cov_beta1, len_CI_beta0 = len_CI_beta0, len_CI_beta1 = len_CI_beta1)
      }, mc.cores = detectCores() - 1) # Use all but one core
      
      # Combine results
      coverages <- do.call(rbind, lapply(simulation_results, function(res) c(res$cov_beta0, res$cov_beta1)))
      length_CI <- do.call(rbind, lapply(simulation_results, function(res) c(res$len_CI_beta0, res$len_CI_beta1)))
      
      beta0_coverage[m, k] <- mean(coverages[, 1])
      beta1_coverage[m, k] <- mean(coverages[, 2])
      
      beta0_CI_length[m, k] <- mean(length_CI[, 1])
      beta1_CI_length[m, k] <- mean(length_CI[, 2])
    }
  }
  list(beta0_coverage = beta0_coverage, beta1_coverage = beta1_coverage, beta0_CI_length = beta0_CI_length, beta1_CI_length = beta1_CI_length)
}

# Run the main function
rt2 <- main_function()
print(rt2)
# Measure execution time using microbenchmark
#timing_results <- microbenchmark(main_function(), times = 10)
#print(timing_results)

