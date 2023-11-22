# Load the required libraries
library(ggplot2)
library(tidyr)
library(LaplacesDemon)


## Private Estimate
get_pi = function(theta, u, w, eps, n){
  set.seed(1)
  w=runif(1, -0.5, 0.5)
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}


# Fudicial Application Setting
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




# CI for true_beta based on B JINI solutions
CI_beta_Fuducial<-function(p_ref, p, B, eps, n_ref, n, alpha, u_ref=runif(n_ref), u=runif(n), w_ref=runif(1, -0.5, 0.5), w=runif(1, -0.5, 0.5), seed=123){
  
  true_beta <- log(p/(1-p)) - log(p_ref/(1-p_ref))
  p_ref_hat<-get_pi(p_ref,u=u_ref, w=w, eps=eps, n=n_ref)
  p_hat<-get_pi(p,u=u, w=w, eps=eps, n=n)
  
  
  #solve the JINI equation B times 
  JINI_solution_p_ref_hat <- JINI_algorithm_Fudicial_Appl(pi0=p_ref_hat, B=B, eps=eps, n=n_ref, seed=seed)
  JINI_solution_p_hat <- JINI_algorithm_Fudicial_Appl(pi0=p_hat, B=B, eps=eps, n=n, seed=seed)
  
  # Calculate the corresponding betas based on the JINI solutions
  beta_hat <-log(JINI_solution_p_hat/(1-JINI_solution_p_hat)) - log(JINI_solution_p_ref_hat/(1-JINI_solution_p_ref_hat))
  
  #calculate the CI for true_betas using the distribution of betas obtained above
  CI_JINI_One_sample_proportion  <- c(quantile(beta_hat, alpha/2, na.rm = TRUE), quantile(beta_hat, 1-(alpha/2), na.rm=TRUE) )
  
  #return(list(true_beta2,CI_JINI_One_sample_proportion,beta2_hat))
  return(CI_JINI_One_sample_proportion)
  
}

# For each group (other than the ref group) Perform the above H times to obtain H CIs for true_betas; check proportion containing true_betas
Coverage<-function(p_ref, p1, p2, p3, B=1000,H=1000, eps=1, n_ref=30, n1=30, n2=37, n3=30, alpha=0.05, seed=123){
  
  true_beta1 <- log(p1/(1-p1)) - log(p_ref/(1-p_ref))
  true_beta2 <- log(p2/(1-p2)) - log(p_ref/(1-p_ref))
  true_beta3 <- log(p3/(1-p3)) - log(p_ref/(1-p_ref))
  
  
  H_matrix<-matrix(NA,H,6)
  for (h in 1:H) {
    CI1<-CI_beta_Fuducial(p_ref=p_ref, p=p1, B=B, eps=eps/4, n_ref=n_ref, n=n1, alpha=alpha, u_ref=runif(n_ref), u=runif(n1), w_ref=runif(1, -0.5, 0.5), w=runif(1, -0.5, 0.5), seed=seed+ 10*h)
    CI2<-CI_beta_Fuducial(p_ref=p_ref, p=p2, B=B, eps=eps/4, n_ref=n_ref, n=n2, alpha=alpha, u_ref=runif(n_ref), u=runif(n2), w_ref=runif(1, -0.5, 0.5), w=runif(1, -0.5, 0.5), seed=seed+ 10*h)
    CI3<-CI_beta_Fuducial(p_ref=p_ref, p=p3, B=B, eps=eps/4, n_ref=n_ref, n=n3, alpha=alpha, u_ref=runif(n_ref), u=runif(n3), w_ref=runif(1, -0.5, 0.5), w=runif(1, -0.5, 0.5), seed=seed+ 10*h)
    
    H_matrix[h,1]<- CI1[1] #lower limit of CI1
    H_matrix[h,2]<- CI1[2]  #upper limit of CI1
    H_matrix[h,3]<- CI2[1] #lower limit of CI2
    H_matrix[h,4]<- CI2[2]  #upper limit of CI2
    H_matrix[h,5]<- CI3[1] #lower limit of CI3
    H_matrix[h,6]<- CI3[2]  #upper limit of CI3
  }
  
  # coverage for each group
  proportion_CI_containing_true_beta1 <- mean(H_matrix[, 1] <= true_beta1 & H_matrix[, 2] >= true_beta2)
  proportion_CI_containing_true_beta2 <- mean(H_matrix[, 3] <= true_beta2 & H_matrix[, 4] >= true_beta2)
  proportion_CI_containing_true_beta3 <- mean(H_matrix[, 5] <= true_beta3 & H_matrix[, 6] >= true_beta2)
  
  
  # Output the proportions for each group
  proportions <- c(Group1 = proportion_CI_containing_true_beta1, 
                   Group2 = proportion_CI_containing_true_beta2, 
                   Group3 = proportion_CI_containing_true_beta3)
  
  return(proportions)
  
  #return(list(true_beta1, true_beta2, true_beta3, proportion_CI_containing_true_beta1,proportion_CI_containing_true_beta2, proportion_CI_containing_true_beta3))
  #return(list(proportion_CI_containing_true_beta1,proportion_CI_containing_true_beta2, proportion_CI_containing_true_beta3))
}




# Ploting coverage for each for different sets of p_ref,p1,p2,p3



# Assume p_values_sets is a list of lists, where each inner list contains p1, p2, p3, p4
plot_CI_betas <- function(p_values_sets, B, H, eps, n_ref, n1, n2, n3, alpha, seed) {
  proportions_data <- data.frame(Group = character(), Proportion = numeric(), Set = character())
  
  for (i in seq_along(p_values_sets)) {
    set <- p_values_sets[[i]]
    p_ref <- set[1]
    p1 <- set[2]
    p2 <- set[3]
    p3 <- set[4]
    
    # Format the set of p-values for labeling
    set_label <- paste("Probabilities:", paste(set, collapse = ", "))
    
    # Call CI_betas function for each set
    proportions <- Coverage(p_ref = p_ref, p1 = p1, p2 = p2, p3 = p3, B = B, H = H, eps = eps, n_ref = n_ref, n1 = n1, n2 = n2, n3 = n3, alpha = alpha, seed = seed)
    
    # Store the results with the formatted set label
    temp_data <- data.frame(Group = names(proportions), Proportion = unlist(proportions), Set = set_label)
    proportions_data <- rbind(proportions_data, temp_data)
  }
  
  # Plot the proportions for each group
  ggplot(proportions_data, aes(x = Set, y = Proportion, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "orange") +
    theme_minimal() +
    labs(title = "Coverage for log odd ratios with Gr1 as ref (H=B=10^3, eps=1,n1=20,n2=15,n3=17,n4=31)", x = "Set of probabilities", y = "Coverage") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Correct way to rotate x-axis labels
}






# Example usage with different sets of values of p 
p_values_sets <- list(
  c(0.3, 0.6, 0.8, 0.4),
  c(0.6, 0.95, 0.45, 0.65),
  c(0.95, 0.2, 0.9, 0.5)
)

plot_CI_betas(p_values_sets, 1000, 1000, 1, 20, 15, 17, 31, 0.05, 1234)






























# 
# 
# plot_CI_betas <- function(p_values_sets, B, H, eps, n_ref, n1, n2, n3, alpha, seed) {
#   proportions_data <- data.frame(Group = character(), Proportion = numeric(), SetIndex = numeric(), SetLabel = character())
# 
#   for (i in seq_along(p_values_sets)) {
#     set <- p_values_sets[[i]]
#     p_ref <- set[1]
#     p1 <- set[2]
#     p2 <- set[3]
#     p3 <- set[4]
# 
#     # Call Coverage function for each set
#     proportions <- Coverage(p_ref = p_ref, p1 = p1, p2 = p2, p3 = p3, B = B, H = H, eps = eps, n_ref = n_ref, n1 = n1, n2 = n2, n3 = n3, alpha, seed)
# 
#     # Store the results
#     temp_data <- data.frame(Group = names(proportions), Proportion = unlist(proportions), SetIndex = i, SetLabel = paste("Set", i))
#     proportions_data <- rbind(proportions_data, temp_data)
#   }
# 
#   # Plot the proportions for each group using line plot
#   ggplot(proportions_data, aes(x = SetIndex, y = Proportion, group = Group, color = Group)) +
#     geom_line() +
#     geom_point() +
#     geom_hline(yintercept = 0.95, linetype = "dashed", color = "orange") +
#     scale_x_continuous(breaks = 1:length(p_values_sets), labels = sapply(p_values_sets, function(x) paste0("Set (", x[1], ", ", x[2], ", ", x[3], ", ", x[4], ")"))) +
#     theme_minimal() +
#     labs(title = "Coverage for Betas", x = "Set of Beta Values", y = "Coverage") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
# }
# 
# 
# # Example usage with different sets of values of p 
# p_values_sets <- list(
#   c(0.3, 0.6, 0.8, 0.4),
#   c(0.6, 0.95, 0.45, 0.65),
#   c(0.95, 0.2, 0.9, 0.5)
# )
# 
# plot_CI_betas(p_values_sets, 1000, 1000, 1, 20, 15, 17, 31, 0.05, 1234)
# 
# 
# 
