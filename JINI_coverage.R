
# Load the required libraries
library(ggplot2)
library(tidyr)


# Define a function for the interpolation
inter <- function(data, p) {
  sorted_data <- sort(data)  # Sort the data in ascending order
  n <- length(sorted_data)   # Number of data points
  
  # Find the smallest rank k such that F(U_(k)) >= p. Note: F(U_(k))=k/n
  k <- 1
  while (k <= n && k / n < p) {
    k <- k + 1
  }
  
  # If p falls between two ranks, perform linear interpolation
  if (k < n && k > 1 && k / n > p) {
    u_k_minus_1 <- sorted_data[k - 1]
    u_k <- sorted_data[k]
    u_k_plus_1 <- sorted_data[k+1]
    v_k <- (u_k + u_k_plus_1)/4
    v_k_minus <- (u_k_minus_1 + u_k)/4
    F_u_k_minus_1 <- (k - 1) / n
    F_u_k <- k / n
    
    inter_value <- u_k - ((u_k - u_k_minus_1)* (F_u_k-p)) / ((F_u_k - F_u_k_minus_1))
  } else {
    
    if (k <= n) {
      inter_value <- sorted_data[k]
    } 
    
  }
  
  return(inter_value)
}







## Private Estimate
get_pi = function(theta, u, w, eps, n){
  set.seed(1)
  w=runif(1, -0.5, 0.5)
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}


# JINI Algorithm (Via Interpolation)

JINI_Inter <- function(pi0, B, eps, n, seed){
  JINI_solution <- numeric(B) 
  
  
  for (i in 1:B) {
    set.seed(seed*i)
    Wj <- runif(1, -0.5, 0.5)
    simulated_sample = runif(n, 0, 1)
    pi_j_star <- pi0 - 1/(eps*n)*sign(Wj)*log(1 - 2*abs(Wj))
    if (pi_j_star >= 1) {
      JINI_solution[i] <- 1#max(simulated_sample)
    } else if (pi_j_star <= 0) {
      JINI_solution[i] <- 0#min(simulated_sample)
    } else {
      JINI_solution[i] <- inter(simulated_sample, pi_j_star)
    }
  }
  JINI_solution
}




# Fudicial Setting
JINI_algorithm_Fudicial_Appl <-function(pi0, B, eps, n, seed){
  res = rep(NA, B)
  counter = 0
  
  while (counter < B){
    set.seed(seed*counter)
    Wj = runif(1, -0.5, 0.5)
    pi_j_star = pi0 - 1/(eps*n)*sign(Wj)*log(1 - 2*abs(Wj))
    
    if (pi_j_star < 1 & pi_j_star > 0){
      counter = counter + 1
      thetaj = rbeta(1, n*pi_j_star+0.5, n*(1-pi_j_star)+0.5, ncp = 0)
      res[counter] = thetaj
    }else if(pi_j_star >= 1){
      counter = counter + 1
      res[counter] = 1
    }else {
      counter = counter + 1
      res[counter] = 0
    }
  }
  res
}







# CI for One-sample Test of Proportion Via Interpolation (Application Setting)
One_sample_CI_application<-function(theta0=0.7, B=100, eps=1, n=30, alpha=0.05, u=runif(n), w=runif(1, -0.5, 0.5), seed=123){
  pi0<-get_pi(theta0,u=u, w=w, eps=eps, n=n)
  #solve the JINI equation B times save the solution in JINI_solution
  JINI_solution <- JINI_Inter(pi0=pi0, B=B, eps=eps, n=n, seed=seed)
  
  #calculate the CI using the solution obtained above
  CI_JINI_One_sample_proportion  <- c(quantile(JINI_solution,   alpha/2), quantile(JINI_solution, 1-(alpha/2)))
  
  return(CI_JINI_One_sample_proportion)
  
}





# CI for One-sample Test of Proportion Via Fudicial (Application Setting)
One_sample_CI_Fiducial_Appl<-function(theta0=0.7, B=100, eps=1, n=30, alpha=0.05, u=runif(n), w=runif(1, -0.5, 0.5), seed=123){
  pi0<-get_pi(theta0,u=u, w=w, eps=eps, n=n)
  
  #solve the JINI equation B times save the solution in JINI_solution
  JINI_solution <- JINI_algorithm_Fudicial_Appl(pi0=pi0, B=B, eps=eps, n=n, seed=seed)
  
  #calculate the CI using the solution obtained above
  CI_JINI_One_sample_proportion  <- c(quantile(JINI_solution,   alpha/2), quantile(JINI_solution, 1-(alpha/2)))
  
  return(CI_JINI_One_sample_proportion)
  
}







Non_private_CI<-function(p,n,alpha){
  # Simulate Bernoulli trials
  observed_sample<- rbinom(n, 1, p)
  success_proportion <- sum(observed_sample) / n
  
  # Compute 95% confidence interval (using normal approximation)
  se <- sqrt(success_proportion * (1 - success_proportion) / n)
  margin_of_error <- qnorm(1 - alpha / 2) * se
  ci_lower <- success_proportion - margin_of_error
  ci_upper <- success_proportion + margin_of_error
  
  return(c(ci_lower, ci_upper))
  
  
}




# Coverage for One-Sample Test of Proportion






CI_coverage_comparison <- function(B = 100, thetavalues = seq(from = 0, to = 1, by = 0.1),
                                   eps = 1, n = 30, seed = 123, H = 1000, alpha = 0.05) {
  
  NonPriv_Lower <- numeric(H)
  NonPriv_Upper <- numeric(H)
  Priv_Lower <- numeric(H)
  Priv_Upper <- numeric(H)
  Priv_Lower_Fudicial <- numeric(H)
  Priv_Upper_Fudicial <- numeric(H)
  
  
  
  # Helper function to compute the coverage probabilities for a given true proportion (p)
  CI_proportion_covered <- function(p) {
    
    for (h in 1:H) {
      set.seed(h*10)
      
      # Nonprivate
      CI_NP<-Non_private_CI(p=p,n=n,alpha = alpha)
      NonPriv_Lower[h] <- CI_NP[1]
      NonPriv_Upper[h] <- CI_NP[2]
      
      # Private Lower and Upper CI Via Interpolation and Fudicial
      Priv_Lower[h] <- One_sample_CI_application(B = B, alpha = alpha, n = n, theta0 = p, eps = eps, u=runif(n), w=runif(1, -0.5, 0.5), seed = seed + h)[1]
      Priv_Upper[h] <- One_sample_CI_application(B = B, alpha = alpha, n = n, theta0 = p, eps = eps, u=runif(n), w=runif(1, -0.5, 0.5), seed = seed + h)[2]
      Priv_Lower_Fudicial[h] <- One_sample_CI_Fiducial_Appl(B = B, alpha = alpha, n = n, theta0 = p, eps = eps, u=runif(n), w=runif(1, -0.5, 0.5), seed = seed + h)[1]
      Priv_Upper_Fudicial[h] <- One_sample_CI_Fiducial_Appl(B = B, alpha = alpha, n = n, theta0 = p, eps = eps,u=runif(n), w=runif(1, -0.5, 0.5), seed = seed + h)[2]
    }
    
    # Checking proportion of the H CIs that contains p
    NonPriv_coverage_prop <- mean(NonPriv_Lower < p & NonPriv_Upper > p)
    Priv_coverage_prop <- mean(Priv_Lower < p & Priv_Upper > p)
    Priv_coverage_prop_Fudicial <- mean(Priv_Lower_Fudicial < p & Priv_Upper_Fudicial > p)
    
    return(c(NonPriv_coverage_prop, Priv_coverage_prop, Priv_coverage_prop_Fudicial))
  }
  
  # Computing coverage probabilities for different values of thetavalues
  coverage_probs_NonPriv <- numeric(length(thetavalues))
  coverage_probs_Priv <- numeric(length(thetavalues))
  coverage_probs_Priv_Fudicial <- numeric(length(thetavalues))
  
  for (j in seq_along(thetavalues)) {
    coverage_probs_NonPriv[j] <- CI_proportion_covered(thetavalues[j])[1]
    coverage_probs_Priv[j] <- CI_proportion_covered(thetavalues[j])[2]
    coverage_probs_Priv_Fudicial[j] <- CI_proportion_covered(thetavalues[j])[3]
  }
  
  # Creating a data frame to store the coverage probabilities
  df <- data.frame(thetavalues = thetavalues,
                   coverage_NonPriv = coverage_probs_NonPriv,
                   coverage_Priv = coverage_probs_Priv,
                   coverage_Priv_Fudicial = coverage_probs_Priv_Fudicial)
  
  # Converting the data frame from wide to long format for plotting with ggplot2
  df_long <- tidyr::gather(df, key = method, value = coverage, -thetavalues)
  
  # Define the colors and labels for each method
  colors <- c("coverage_NonPriv" = "black", "coverage_Priv" = "red", "coverage_Priv_Fudicial" = "green", "True coverage" = "blue")
  labels <- c("coverage_NonPriv" = "Nonprivate Coverage", "coverage_Priv" = "Private Coverage Interpolation", "coverage_Priv_Fudicial" = "Private Covergae Fudicial", "True coverage" = "True coverage")
  
  # Create the plot with ggplot2
  ggplot(data = df_long, aes(x = thetavalues, y = coverage, color = method)) +
    geom_line() +
    scale_color_manual(values = colors, labels = labels) +
    geom_point() +
    geom_hline(yintercept = 1-alpha, color = "blue") +
    labs(title = "Coverage (B=H=10^3, n=30, eps=1)", x = "thetavalues", y = "Coverage probabilities") +
    ylim(0, 1) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(linetype = 1, shape = 16)))
}


#To see plots of coverage in the window
CI_coverage_comparison(B = 10^3, thetavalues = seq(from = 0, to = 1, by = 0.1), eps = 1, n = 30, seed = 345, H = 10^3, alpha = 0.05)

# # Save the plot as JINI_coverage.png (This is for running on HPC)
# plot_file <- "JINI_coverage.png"
# ggsave(plot = CI_coverage_comparison(B = 10, thetavalues = seq(from = 0, to = 1, by = 0.1), eps = 1, n = 30, seed = 345, iter = 100, alpha = 0.05),
#        filename = plot_file)
