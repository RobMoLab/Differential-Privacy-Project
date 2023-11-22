
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
  if (k <= n && k > 1 && k / n > p) {
    u_k_minus_1 <- sorted_data[k - 1]
    u_k <- sorted_data[k]
    F_u_k_minus_1 <- (k - 1) / n
    F_u_k <- k / n
    
    inter_value <- u_k - ((u_k - u_k_minus_1)* (F_u_k-p)) / ((F_u_k - F_u_k_minus_1))
  } else {
    # If p exactly matches the ECDF value, use the data point directly
    if (k <= n) {
      inter_value <- sorted_data[k]
    } 
    #else {
    # If p exceeds the maximum ECDF value, use the maximum data point
    #  inter_value <- sorted_data[n]
    # }
  }
  
  return(inter_value)
}





## Private Estimate
get_pi = function(theta, u, w, eps, n){
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}


# JINI Algorithm (Via Interpolation)

JINI_Inter <- function(pi0, B, eps, n, seed){
  set.seed(seed)
  JINI_solution <- numeric(B) 
  
  
  for (i in 1:B) {
    Wj <- runif(1, -0.5, 0.5)
    simulated_sample = runif(n, 0, 1)
    pi_j_star <- pi0 - 1/(eps*n)*sign(Wj)*log(1 - 2*abs(Wj))
    if (pi_j_star >= 1) {
      JINI_solution[i] <- 1 #max(simulated_sample)
    } else if (pi_j_star <= 0) {
      JINI_solution[i] <- 0 #min(simulated_sample)
    } else {
      JINI_solution[i] <- inter(simulated_sample, pi_j_star)
    }
  }
  JINI_solution
}






# CI for One-sample Test of Proportion Via Interpolation (Application Setting)

One_sample_CI_application<-function(theta0=0.7, B=100, eps=1, n=30, alpha=0.05, seed=123){
  pi0<-get_pi(theta0,u=runif(n), w=runif(1, -0.5, 0.5), eps=eps, n=n)
  #solve the JINI equation B times save the solution in JINI_solution
  JINI_solution <- JINI_Inter(pi0=pi0, B=B, eps=eps, n=n, seed=seed)
  
  #calculate the CI using the solution obtained above
  CI_JINI_One_sample_proportion  <- c(quantile(JINI_solution,   alpha/2), quantile(JINI_solution, 1-(alpha/2)))
  
  return(CI_JINI_One_sample_proportion)
  
}




# Fudicial Application Setting

JINI_algorithm_Fudicial_Appl <-function(pi0, B, eps, n, seed){
  set.seed(seed)
  res = rep(NA, B)
  counter = 0
  
  while (counter < B){
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

One_sample_CI_application<-function(theta0=0.7, B=100, eps=1, n=30, alpha=0.05, seed=123){
  pi0<-get_pi(theta0,u=runif(n), w=runif(1, -0.5, 0.5), eps=eps, n=n)
  #solve the JINI equation B times save the solution in JINI_solution
  JINI_solution <- JINI_Inter(pi0=pi0, B=B, eps=eps, n=n, seed=seed)
  
  #calculate the CI using the solution obtained above
  CI_JINI_One_sample_proportion  <- c(quantile(JINI_solution,   alpha/2), quantile(JINI_solution, 1-(alpha/2)))
  
  return(CI_JINI_One_sample_proportion)
  
}




# CI for One-sample Test of Proportion Via Fiducial (Application Setting)


One_sample_CI_Fiducial_Appl<-function(theta0=0.7, B=100, eps=1, n=30, alpha=0.05, seed=123){
  pi0<-get_pi(theta0,u=runif(n), w=runif(1, -0.5, 0.5), eps=eps, n=n)
  
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






# Precision of CIs

precision_and_plot <- function(B = 100, theta0 = 0.65, nvalues = c(10, 15, 20, 25, 30, 40, 60, 100, 120, 150, 180, 200, 250, 300, 500),
                               eps = 0.1, n = 30, seed = 123, iter = 1000, alpha = 0.05) {
  precision <- function(n, iter = 100, B = 100, eps = 1, alpha = 0.05, theta0 = 0.65, seed = 123) {
    NonPriv_precision <- numeric(iter)
    Priv_precision <- numeric(iter)
    Priv_precision_Fiducial <- numeric(iter)
    for (i in 1:iter) {
      observed_sample <- runif(n, 0, 1)
      initial_est <- ecdf(observed_sample)(theta0)
      
      # Nonprivate
      #se1 <- sqrt(initial_est * (1 - initial_est) / n)
      #margin_of_error <- qnorm(1 - alpha / 2) * se1
      CI_Np<-Non_private_CI(p=theta0, n=n, alpha = alpha)
      NonPriv_Lower <- CI_Np[1]  
      NonPriv_Upper <-CI_Np[2] 
      NonPriv_precision[i] <- NonPriv_Upper - NonPriv_Lower
      
      # Private
      CI_P_Inter<-One_sample_CI_application(B = B, alpha = alpha, n = n, theta0 = theta0, eps = eps, seed = seed + i)
      CI_P_Fudicial<-One_sample_CI_Fiducial_Appl(B = B, alpha = alpha, n = n, theta0 = theta0, eps = eps, seed = seed + i)
      Priv_Lower <- CI_P_Inter[1]
      Priv_Upper <- CI_P_Inter[2]
      Priv_Lower_Fiducial <- CI_P_Fudicial[1]
      Priv_Upper_Fiducial <- CI_P_Fudicial[2]
      Priv_length <- Priv_Upper - Priv_Lower
      Priv_length_Fudicial <- Priv_Upper_Fiducial - Priv_Lower_Fiducial
      Priv_precision[i] <- Priv_length
      Priv_precision_Fiducial[i] <- Priv_length_Fudicial
    }
    # Average Precision
    Ave_NonPriv_precision <- mean(NonPriv_precision)
    Ave_Priv_precision <- mean(Priv_precision)
    Ave_Priv_precision_Fiducial <- mean(Priv_precision_Fiducial)
    
    return(c(Ave_NonPriv_precision, Ave_Priv_precision, Ave_Priv_precision_Fiducial))
  }
  
  # Calculate precision for each value of n
  precision_list <- lapply(nvalues, precision, iter = iter, B = B, eps = eps, alpha = alpha, theta0 = theta0, seed = seed)
  
  # Separate Ave_NonPriv_precision and Ave_Priv_precision into individual vectors
  Ave_NonPriv_precision <- sapply(precision_list, function(x) x[1])
  Ave_Priv_precision <- sapply(precision_list, function(x) x[2])
  Ave_Priv_precision_Fiducial <- sapply(precision_list, function(x) x[3])
  
  df <- data.frame(nvalues, Ave_NonPriv_precision, Ave_Priv_precision, Ave_Priv_precision_Fiducial)
  
  library(ggplot2)
  plot <- ggplot(df, aes(nvalues)) +
    geom_point(aes(y = Ave_NonPriv_precision, color = "NonPrivate")) +
    geom_line(aes(y = Ave_NonPriv_precision, color = "NonPrivate")) +
    geom_point(aes(y = Ave_Priv_precision, color = "Private")) +
    geom_line(aes(y = Ave_Priv_precision, color = "Private")) +
    geom_point(aes(y = Ave_Priv_precision_Fiducial, color = "Private_Fiducial")) +
    geom_line(aes(y = Ave_Priv_precision_Fiducial, color = "Private_Fiducial")) +
    labs(title = "Length of Private and Non-private CI (theta0=0.2, B=10^3, H=10^4, eps=1) ", x = "Sample size n", y = "Length", color = "Settings") +
    scale_color_manual(values = c("NonPrivate" = "black", "Private Interpolation" = "red", "Private_Fiducial" = "green"))
  
  colnames(df) <- c("sample size n", "NonPrivate CI", "Private CI", "Private CI Fiducial")
  return(list(df = df, plot = plot))
}

# Example usage
result <- precision_and_plot(B = 10^3, theta0  = 0.2, nvalues = c(8*2^(0:6), 700, 900, 1200),
                             eps = 1, seed = 123, iter = 10^4, alpha = 0.05)

#print(result$df) # Display the data frame
print(result$plot) # Display the plot
