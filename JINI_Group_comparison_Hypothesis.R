# Load the required libraries
library(ggplot2)
library(tidyr)
library(LaplacesDemon)


## Private Estimate for four groups
get_pi = function(p1,p2,p3,p4, eps, n1,n2,n3,n4,seed){
  set.seed(seed)
  G1<-mean(rbinom(n1,1,p1)) + rlaplace(1,0,1/(eps*n1))
  G2<-mean(rbinom(n2,1,p2)) + rlaplace(1,0,1/(eps*n2))
  G3<-mean(rbinom(n3,1,p3)) + rlaplace(1,0,1/(eps*n3))
  G4<-mean(rbinom(n4,1,p4)) + rlaplace(1,0,1/(eps*n4))
  return(c(G1,G2,G3,G4))
}



# Fudicial Application 
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

# # # Safe log-odds ratio calculation to avoid division by zero
# safe_log_odds <- function(p_hat, p1_hat) {
#   if (is.na(p_hat) || is.na(p1_hat) || is.nan(p_hat) || is.nan(p1_hat) ||
#       p_hat <= 0 || p_hat >= 1 || p1_hat <= 0 || p1_hat >= 1) {
#     return(NA)
#   } else {
#     return(log(p_hat / (1 - p_hat)) - log(p1_hat / (1 - p1_hat)))
#   }
# }
# 





###################################  Group Comparison       ####################################


# we test the hypotheses; H0: beta >= 0 Vs Ha: beta < 0; 
# where beta is the log odd ratio of Gr2, Gr3 and Gr3 wrt Gr1
#So Group 1 is the reference level


betas_with_testing <- function(p1, p2, p3, p4, B, H, eps, n1, n2, n3, n4, alpha = 0.05, seed) {
  
  # True log odd ratios with Group1 as a reference group
  true_beta2<-log(p2/(1-p2))-log(p1/(1-p1))
  true_beta3<-log(p3/(1-p3))-log(p1/(1-p1))
  true_beta4<-log(p4/(1-p4))-log(p1/(1-p1))
  
  
  # Initialize matrix to store p-values for each group and each iteration
  p_values <-beta_estimates<- matrix(NA, nrow = H, ncol = 3)
  
  for (h in 1:H) {
    set.seed(h)
    
    # Observed privatized proportions for each group
    observed_proportion<-get_pi(p1, p2, p3, p4, eps/4, n1, n2, n3, n4, seed)
    
    observed_p1_hat<-observed_proportion[1]
    observed_p2_hat<-observed_proportion[2]
    observed_p3_hat<-observed_proportion[3]
    observed_p4_hat<-observed_proportion[4]
    
    # Calculate observed betas with Group1 as a reference group
    observed_beta2 <- log(observed_p2_hat/(1-observed_p2_hat)) - log(observed_p1_hat/(1-observed_p1_hat))
    observed_beta3 <- log(observed_p3_hat/(1-observed_p3_hat)) - log(observed_p1_hat/(1-observed_p1_hat))
    observed_beta4 <- log(observed_p4_hat/(1-observed_p4_hat)) - log(observed_p1_hat/(1-observed_p1_hat))
    
    
    # observed_beta2 <- safe_log_odds_ratio(observed_p2_hat, observed_p1_hat)
    # observed_beta3 <- safe_log_odds_ratio(observed_p3_hat, observed_p1_hat)
    # observed_beta4 <- safe_log_odds_ratio(observed_p4_hat, observed_p1_hat)

    
    
    
    
    # Get JINI empirical distribution of p for each group
    emp_dis_p1_hat <- JINI_algorithm_Fudicial_Appl(observed_p1_hat, B, eps/4, n1, seed + 2*h)
    emp_dis_p2_hat <- JINI_algorithm_Fudicial_Appl(observed_p2_hat, B, eps/4, n2, seed + 2*h)
    emp_dis_p3_hat <- JINI_algorithm_Fudicial_Appl(observed_p3_hat, B, eps/4, n3, seed + 2*h)
    emp_dis_p4_hat <-JINI_algorithm_Fudicial_Appl(observed_p4_hat, B, eps/4, n4, seed + 2*h)
    
    
    
    # Calculate empirical distribution of beta for each group
    Jini_emp_dis_beta2 <- log(emp_dis_p2_hat/(1-emp_dis_p2_hat)) - log(emp_dis_p1_hat/(1-emp_dis_p1_hat))
    Jini_emp_dis_beta3 <- log(emp_dis_p3_hat/(1-emp_dis_p3_hat)) - log(emp_dis_p1_hat/(1-emp_dis_p1_hat))
    Jini_emp_dis_beta4 <- log(emp_dis_p4_hat/(1-emp_dis_p4_hat)) - log(emp_dis_p1_hat/(1-emp_dis_p1_hat))
    
    # Jini_emp_dis_beta2 <- sapply(emp_dis_p2_hat, safe_log_odds, p1_hat = emp_dis_p1_hat)
    # Jini_emp_dis_beta3 <- sapply(emp_dis_p3_hat, safe_log_odds, p1_hat = emp_dis_p1_hat)
    # Jini_emp_dis_beta4 <- sapply(emp_dis_p4_hat, safe_log_odds, p1_hat = emp_dis_p1_hat)

    
    

    
    # Calculate p-values for this iteration
    p_values[h, 1] <- sum(Jini_emp_dis_beta2 >=0, na.rm = TRUE)/(B+1) 
    p_values[h, 2] <- sum(Jini_emp_dis_beta3 >=0, na.rm = TRUE)/(B+1) 
    p_values[h, 3] <- sum(Jini_emp_dis_beta4 >=0, na.rm = TRUE)/(B+1)

    
    
    
    
    # Calculate average of B estimated betas for this iteration
    beta_estimates[h, 1] <- mean(Jini_emp_dis_beta2, na.rm = TRUE) 
    beta_estimates[h, 2] <- mean(Jini_emp_dis_beta3, na.rm = TRUE) 
    beta_estimates[h, 3] <- mean(Jini_emp_dis_beta4, na.rm = TRUE) 
    
    
    
    
  }
  
  #beta_estimates[is.infinite(beta_estimates)] <- NA
 
  # Handle Inf and -Inf values
  beta_estimates[is.infinite(beta_estimates)] <- NA
  beta_estimates[is.nan(beta_estimates)] <- NA
  
  
  # Aggregate results and calculate mean p-values for each group
  
  test_results <- data.frame(
    group = c("Group2", "Group3", "Group4"),
    #beta_estimates_mean = colMeans(beta_estimates,)
    beta_estimates_mean = colMeans(beta_estimates, na.rm = TRUE),
    true_beta = c(true_beta2,true_beta3,true_beta4),
    p_value_mean = colMeans(p_values, na.rm = TRUE),
    power = colMeans(p_values<0.05, na.rm = TRUE)
  )
  
  # Add decision column based on mean p-values
  test_results$reject_H0 <- test_results$p_value_mean < alpha
  return(list(test_results = test_results))
  
}

# Example: Comparing four groups
p1 <- 0.7  # Example probabilities for each group
p2 <- 0.8
p3 <- 0.65     
p4 <- 0.35
B <- 10^3  # Number of bootstrap samples
H <- 10^3  # Number of iterations
eps <- 1  
n1 <- 200
n2 <- 100
n3 <- 120
n4 <- 100  
alpha <- 0.05 
seed <- 12345

# Run the hypothesis testing
test_results <- betas_with_testing(p1, p2, p3, p4, B, H, eps, n1, n2, n3, n4, alpha, seed)
print(test_results)























