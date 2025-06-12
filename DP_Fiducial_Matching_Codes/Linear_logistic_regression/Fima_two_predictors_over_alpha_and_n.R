# Load necessary package
library(microbenchmark)

# Function Definitions
inverse_logit <- function(nu) {
  exp(nu) / (1 + exp(nu))
}


# private proportion
get_pi <- function(p_hat, eps, n) {
  w = runif(1, -0.5, 0.5)
  # Ensure w is strictly between -0.5 and 0.5
  if (w == -0.5) w <- -0.5 + .Machine$double.eps
  if (w == 0.5) w <- 0.5 - .Machine$double.eps
  p_hat - 1 / (eps * n) * sign(w) * log(1 - 2 * abs(w))
}

# JINI Fiducial Algorithm for one proportion
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

# Encapsulate main code in a function
main_function <- function() {
  # Parameters
  alpha_values <- c(0.01, 0.05, 0.1)
  N <- c(15, 30, 100, 200, 500, 1000, 2000)
  eps <- 1
  H <- 10^4
  beta0 <- 0.5
  beta1 <- 2
  beta2 <- -2
  beta <- c(beta0, beta1, beta2)
  p1 <- 0.5
  p2 <- 0.5
  seed <- 144
  # Results Storage
  initialize_matrices <- function(dim1, dim2) {
    matrix(NA, dim1, dim2)
  }
  
  beta0_coverage <- initialize_matrices(length(N), length(alpha_values))
  beta1_coverage <- initialize_matrices(length(N), length(alpha_values))
  beta2_coverage <- initialize_matrices(length(N), length(alpha_values))
  beta0_CI_length <- initialize_matrices(length(N), length(alpha_values))
  beta1_CI_length <- initialize_matrices(length(N), length(alpha_values))
  beta2_CI_length <- initialize_matrices(length(N), length(alpha_values))
  
  # Main Loops
  for (k in 1:length(alpha_values)) {
    for (m in 1:length(N)) {
      set.seed(m + seed)
      X1 <- rbinom(N[m], 1, p1)
      X2 <- rbinom(N[m], 1, p2)
      nu <- beta0 + beta1 * X1 + beta2 * X2
      probs <- inverse_logit(nu)
      seed <- N[m] + seed 
      
      # Approximation
      W1 = (X1 == 0) *(X2 == 0)
      W2 = (X1 == 1) *(X2 == 0)
      W3 = (X1 == 0)*(X2 == 1)
      
      # Simulation
      B <- 10^4
      coverages <- matrix(NA, B, length(beta))
      length_CI <- matrix(NA, B, length(beta))
      
      for (i in 1:B) {
        set.seed(i)
        u <- runif(N[m])
        Y <- as.numeric(probs > u)
        
        pi1 = mean(Y[W1==1])
        pi2 = mean(Y[W2==1])
        pi3 = mean(Y[W3==1])
        
        n1 = sum(W1)
        n2 = sum(W2)
        n3 = sum(W3)
        
        pi1_priv <- get_pi(pi1, eps, n1)
        pi2_priv <- get_pi(pi2, eps, n2)
        pi3_priv <- get_pi(pi3, eps, n3)
        
        jini_dis_pi1 <- JINI_algorithm_Fiducial(pi1_priv, H, eps, n1, seed=seed + 2 * i)
        jini_dis_pi2 <- JINI_algorithm_Fiducial(pi2_priv, H, eps, n2, seed=seed + 2 * i)
        jini_dis_pi3 <- JINI_algorithm_Fiducial(pi3_priv, H, eps, n3, seed=seed + 2 * i)
        
        beta_0 <- log(jini_dis_pi1 / (1 - jini_dis_pi1))
        beta_1 <- log(jini_dis_pi2 / (1 - jini_dis_pi2)) - beta_0
        beta_2 <- log(jini_dis_pi3 / (1 - jini_dis_pi3)) - beta_0
        
        CI_beta0 <- quantile(beta_0, c(alpha_values[k] / 2, 1 - alpha_values[k] / 2))
        CI_beta1 <- quantile(beta_1, c(alpha_values[k] / 2, 1 - alpha_values[k] / 2))
        CI_beta2 <- quantile(beta_2, c(alpha_values[k] / 2, 1 - alpha_values[k] / 2))
        
        len_CI_beta0 <- CI_beta0[2] - CI_beta0[1]
        len_CI_beta1 <- CI_beta1[2] - CI_beta1[1]
        len_CI_beta2 <- CI_beta2[2] - CI_beta2[1]
        
        cov_beta0 <- (CI_beta0[1] < beta0) & (CI_beta0[2] > beta0)
        cov_beta1 <- (CI_beta1[1] < beta1) & (CI_beta1[2] > beta1)
        cov_beta2 <- (CI_beta2[1] < beta2) & (CI_beta2[2] > beta2)
        
        coverages[i, ] <- c(cov_beta0, cov_beta1, cov_beta2)
        length_CI[i, ] <- c(len_CI_beta0, len_CI_beta1, len_CI_beta2)
      }
      
      beta0_coverage[m, k] <- mean(coverages[, 1])
      beta1_coverage[m, k] <- mean(coverages[, 2])
      beta2_coverage[m, k] <- mean(coverages[, 3])
      
      beta0_CI_length[m, k] <- mean(length_CI[, 1])
      beta1_CI_length[m, k] <- mean(length_CI[, 2])
      beta2_CI_length[m, k] <- mean(length_CI[, 3])
    }
  }
  
  list(beta0_coverage = beta0_coverage, beta1_coverage = beta1_coverage, beta2_coverage = beta2_coverage, beta0_CI_length = beta0_CI_length, beta1_CI_length = beta1_CI_length, beta2_CI_length = beta2_CI_length)
}

# Run the main function
rt4 <- main_function()
print(rt4)



# # Measure execution time using microbenchmark
# timing_results <- microbenchmark(main_function(), times = 5)
# print(timing_results)

# Save timing_results as .rda file
#save(timing_results, file = file.path(output_dir, "timing_results_fiducial_two_predictor_alpha_n.rda"))
