# Privatize Proportions
get_pi = function(theta, u, w, eps, n){
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}


##### Fiducial Approach to Solve JINI  #########

# Private
JINI_algorithm_Fiducial_Appl <-function(pi0, B, eps, n, seed){
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


# JINI Distribution of beta
jini_logistic_fiducial_Priv = function(pi1, pi2, B, eps, n1, n2, seed){
  set.seed(seed)
  res = rep(NA, B)
  counter = 0
  
  while (counter < B){
    distri_grp1 = JINI_algorithm_Fiducial_Appl(pi0 = pi1, B=B, eps = eps/2, n = n1, seed=seed)
    distri_grp2 = JINI_algorithm_Fiducial_Appl(pi0 = pi2, B=B, eps = eps/2, n = n2, seed=seed)
    
    counter = counter + 1
    res[counter] = log(distri_grp2/(1-distri_grp2)) - log(distri_grp1/(1-distri_grp1)) 
  }
  res
}



# Record the start time
#start_time <- Sys.time()


############ Level Analysis  ####################
#      beta = log(p2/(1-p2)) - log(p1/(1-p1))

# H0: beta = 0      Vs    H1: beta < 0

H = 1
B = 1000
n1 = 30 
n2 = 30 
theta0 = c(0.1, 0.3, 0.5, 0.7, 0.9)
eps = 2
res_PJ2 = matrix(NA, H, length(theta0))

for (j in 1:length(theta0)){
  for (i in 1:H){
    set.seed(i)
    
    # Privatized proportions
    pi1 = get_pi(theta = theta0[j], u = runif(n1), w = runif(1, -0.5, 0.5), eps = eps/2, n = n1)
    pi2 = get_pi(theta = theta0[j], u = runif(n2), w = runif(1, -0.5, 0.5), eps = eps/2, n = n2)
    
    
    # Distribution of beta
    emp_density_PJ = jini_logistic_fiducial_Priv(pi1 = pi1, pi2 = pi2, B = B, eps = eps, n1 = n1, n2 = n2, seed = i + 2*B)
    
    
    # P-values
    res_PJ2[i,j] = (sum(emp_density_PJ < 0)+1)/(B + 1)
    
  }
}


level_res_PJ2 <- apply(res_PJ2 < 0.05, 2, mean, na.rm = TRUE)

# #To run on computer with output
# plot(theta0, level_res_PJ2 , type = "b", pch = 19, col = "blue", xlab = "theta values", ylab = "levels", ylim = c(0, 1), main = "Level of Private Logistic Regression (H=B=10^3, eps=1,alpha=0.05, n1=30, n2=40)")
# abline(h = 0.05, col = "purple")


#To run on HPC
# Plot of level
pdf("level_Private_Logistic_JINI.pdf")  # Save plot as PDF on Auburn HPC
plot(theta0, level_res_PJ2 , type = "b", pch = 19, col = "blue", xlab = "theta values", ylab = "levels", ylim = c(0, 1), main = "Level of Two Samples Private JINI Test 
     (H=B=10^3, eps=1,alpha=0.05, n1=30, n2=40)")
abline(h = 0.05, col = "purple")
dev.off()




# # Record the end time
# end_time <- Sys.time()
# 
# # Calculate the execution time
# execution_time <- end_time - start_time
# print(paste("Execution time:", execution_time))

