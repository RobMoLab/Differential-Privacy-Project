# Load necessary package
library(microbenchmark)



    ################# Helper functions ################
# inverse logistic function
inverse_logit <- function(nu) {
  exp(nu) / (1 + exp(nu))
}


# Function for Privatizing proportions
get_pi <- function(p_hat, eps, n) {
  w = runif(1, -0.5, 0.5)
  # Ensure w is strictly between -0.5 and 0.5
  if (w == -0.5) w <- -0.5 + .Machine$double.eps
  if (w == 0.5) w <- 0.5 - .Machine$double.eps
  p_hat - 1 / (eps * n) * sign(w) * log(1 - 2 * abs(w))
}


# Fima Algorithm for one proportion
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




# Simulation for one predictor



# Encapsulate main code in a function for execution time computation
main_function <- function() {
  
  # Parameters
  alpha_values <- c(0.01, 0.05, 0.1)          # significant levels
  N <- c(15, 30, 100, 200, 500, 1000, 2000)   # sample sizes
  eps <- 1
  H <- 10^4
  beta0 <- 0.5
  beta1 <- 2
  beta <- c(beta0, beta1)
  p1 <- 0.5                                     # parameter for generating X1 from Bern(p1)
  
  # Results Storage
  initialize_matrices <- function(dim1, dim2) {
    matrix(NA, dim1, dim2)
  }
  
  beta0_coverage <- initialize_matrices(length(N), length(alpha_values))
  beta1_coverage <- initialize_matrices(length(N), length(alpha_values))     #storage for beta1 coverage
  beta0_CI_length <- initialize_matrices(length(N), length(alpha_values))
  beta1_CI_length <- initialize_matrices(length(N), length(alpha_values))   #storage for length of CI for beta1
  
  # Main Loops
  for (k in 1:length(alpha_values)) {
    for (m in 1:length(N)) {
      set.seed(m + 112)
      X1 <- rbinom(N[m], 1, p1)              # for each sample size, X1 is fixed
      nu <- beta0 + beta1 * X1               # model; logit(P(Y=1))
      probs <- inverse_logit(nu)            # P(Y=1)
      seed <- N[m] + 112 
      
      # Simulation
      B <- 10^4
      coverages <- matrix(NA, B, length(beta))
      length_CI <- matrix(NA, B, length(beta))
      
      for (i in 1:B) {
        set.seed(i)
        u <- runif(N[m])            # seed for sampling Y from the model
        Y <- as.numeric(probs > u)  # Sampling Y from the model in the ith replication
        
        # Private JINI Fiducial
        pi0 <- mean(Y[X1 == 0])      # this is pi0 = #{Y=1,X1=0}/#{X1=0}
        pi1 <- mean(Y[X1 == 1])      # this is pi1 = #{Y=1,X1=1}/#{X1=1}
        
        n1 <- sum(X1 == 0)
        n2 <- sum(X1 == 1)
        
        pi0_priv <- get_pi(pi0, eps, n1)    # this is pi0_hat=pi0+Y0
        pi1_priv <- get_pi(pi1, eps, n2)     # this is pi1_hat = pi1+Y1
        
        jini_dis_pi0 <- JINI_algorithm_Fiducial(pi0_priv, H, eps, n1, seed=seed + 2 * i)  # Solving JINI H times 
        jini_dis_pi1 <- JINI_algorithm_Fiducial(pi1_priv, H, eps, n2, seed=seed + 2 * i)
        
        beta_0 <- log(jini_dis_pi0/ (1 - jini_dis_pi0))                           # Distribution of beta0
        beta_1 <- log(jini_dis_pi1 / (1 - jini_dis_pi1)) - beta_0                  # Distribution of beta1  
        
        CI_beta0 <- quantile(beta_0, c(alpha_values[k] / 2, 1 - alpha_values[k] / 2))   # Coverage for beta0 for the ith iteration
        CI_beta1 <- quantile(beta_1, c(alpha_values[k] / 2, 1 - alpha_values[k] / 2))   # Coverage for beta1 for the ith iteration
        
        len_CI_beta0 <- CI_beta0[2] - CI_beta0[1]                           # length of CI for beta0 for the ith iteration
        len_CI_beta1 <- CI_beta1[2] - CI_beta1[1]                           # length of CI for beta1 for the ith iteration
        
        cov_beta0 <- (CI_beta0[1] < beta0) & (CI_beta0[2] > beta0)
        cov_beta1 <- (CI_beta1[1] < beta1) & (CI_beta1[2] > beta1)
        
        coverages[i, ] <- c(cov_beta0, cov_beta1)
        length_CI[i, ] <- c(len_CI_beta0, len_CI_beta1)
      }
      
      beta0_coverage[m, k] <- mean(coverages[, 1])                         # Average coverage for beta0 over B coverages
      beta1_coverage[m, k] <- mean(coverages[, 2])                          # Average coverage for beta1 over B coverages
      
      beta0_CI_length[m, k] <- mean(length_CI[, 1])                       # Average length of CI for beta0 out of B lengths
      beta1_CI_length[m, k] <- mean(length_CI[, 2])                       # Average length of CI for beta1 out of B lengths
    }
  }
  list(beta0_coverage = beta0_coverage, beta1_coverage = beta1_coverage, beta0_CI_length = beta0_CI_length, beta1_CI_length = beta1_CI_length)
}

# Run the main function and store the results
rt1 <- main_function()
print(rt1)


