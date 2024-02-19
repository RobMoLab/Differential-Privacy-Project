get_pi = function(theta, u, w, eps, n){
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}




#### Normal Approximation Approache to Solve JINI
get_theta_tilde = function(pi0, Zi, Wj, eps, n){
  delta = pi0 - 1/(eps*n)*sign(Wj)*log(1 - 2*abs(Wj))
  inter = Zi^2 - 4*n*delta^2 + 4*n*delta
  if (inter < 0){
    sol1 = 0
    sol2 = 1
    objs = abs(c((sqrt((sol1*(1 - sol1))/n)*Zi + sol1 - delta)^2,
                 (sqrt((sol2*(1 - sol2))/n)*Zi + sol2 - delta)^2))
    return(c(0,1)[which.min(objs)])
  }else{
    sol1 = (2*delta*n + Zi^2 + Zi*(Zi^2 - 4*n*delta^2 + 4*n*delta)^(1/2))/(2*(Zi^2 + n))
    sol2 = (2*delta*n + Zi^2 - Zi*(Zi^2 - 4*n*delta^2 + 4*n*delta)^(1/2))/(2*(Zi^2 + n))
    
    objs = abs(c((sqrt((sol1*(1 - sol1))/n)*Zi + sol1 - delta)^2,
                 (sqrt((sol2*(1 - sol2))/n)*Zi + sol2 - delta)^2))
    
    sol = c(sol1, sol2)[which.min(objs)]
    
    if (sol < 0 && sol > 1){
      sol1 = 0
      sol2 = 1
      objs = abs(c((sqrt((sol1*(1 - sol1))/n)*Zi + sol1 - delta)^2,
                   (sqrt((sol2*(1 - sol2))/n)*Zi + sol2 - delta)^2))
      return(c(0,1)[which.min(objs)])
    }else{
      return(sol)
    }
  }
}

jini2samples = function(pi1, pi2, B, eps, n1, n2, seed){
  set.seed(seed)
  res = rep(NA, B)
  counter = 0
  
  while (counter < B){
    inter1 = get_theta_tilde(pi0 = pi1, Zi = rnorm(1), Wj = runif(1, -0.5, 0.5), eps = eps, n = n1)
    inter2 = get_theta_tilde(pi0 = pi2, Zi = rnorm(1), Wj = runif(1, -0.5, 0.5), eps = eps, n = n2)
    
    counter = counter + 1
    res[counter] = inter1 - inter2
  }
  res
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


jini2samples_fiducial_Priv = function(pi1, pi2, B, eps, n1, n2, seed){
  set.seed(seed)
  res = rep(NA, B)
  counter = 0
  
  while (counter < B){
    distri_grp1 = JINI_algorithm_Fiducial_Appl(pi0 = pi1, B=B, eps = eps, n = n1, seed=seed)
    distri_grp2 = JINI_algorithm_Fiducial_Appl(pi0 = pi2, B=B, eps = eps, n = n2, seed=seed)
    
    counter = counter + 1
    res[counter] = distri_grp1 - distri_grp2
  }
  res
}

# Non-private JINI
get_pi_NP = function(theta,  u){
  mean(u < theta) 
}


JINI_NP<-function(pi0,B,n){
  Np_JINI_Fudicial_solution<-numeric(B)
  for (i in 1:B) {
    Np_JINI_Fudicial_solution[i]<- rbeta(1, n * pi0 + 0.5, n * (1 - pi0) + 0.5)
  }
  return(Np_JINI_Fudicial_solution)
}


jini2samples_fiducial_NP = function(pi1, pi2, B, n1, n2, seed){
  set.seed(seed)
  res = rep(NA, B)
  counter = 0
  
  while (counter < B){
    distri_grp1 = JINI_NP(pi0 = pi1, B=B, n = n1)
    distri_grp2 = JINI_NP(pi0 = pi2, B=B, n = n2)
    
    counter = counter + 1
    res[counter] = distri_grp1 - distri_grp2
  }
  res
}

# # Record the end time
# start_time <- Sys.time()

############ Power


H = 10^3
B = 10^3 
n1 = 16*2^(0:5)
n2 = 12*2^(0:5)
theta1 = 0.8
theta2 = 0.9
eps = 1
res_Normal_app2 = res_NPJ2=res_PJ2=res_np2=matrix(NA, H, length(n1))
#res_np2 = matrix(NA, H, length(n1))

for (j in 1:length(n1)){
  
  for (i in 1:H){
    
    set.seed(i)
    pi1 = get_pi(theta = theta1, u = runif(n1[j]), w = runif(1, -0.5, 0.5), eps = eps, n = n1[j])
    pi2 = get_pi(theta = theta2, u = runif(n2[j]), w = runif(1, -0.5, 0.5), eps = eps, n = n2[j])
    pi1_NPJ = get_pi_NP(theta = theta1, u = runif(n1[j]))
    pi2_NPJ = get_pi_NP(theta = theta2, u = runif(n2[j]))
    
    
    #Solving JINI B times 
    # Distribution of JINI with normal approximation
    emp_density_Norm_ap = jini2samples(pi1 = pi1, pi2 = pi2, B = B, eps = eps, n1 = n1[j], n2 = n2[j], seed = i + 2*B)
    
    # Distribution of Private JINI Fiducial approach
    emp_density_PJ = jini2samples_fiducial_Priv(pi1 = pi1, pi2 = pi2, B = B, eps = eps, n1 = n1[j], n2 = n2[j], seed = i + 2*B)
    
    # Distribution of Non-private JINI Fiducial approach
    emp_density_NPJ = jini2samples_fiducial_NP(pi1 = pi1_NPJ, pi2 = pi2_NPJ, B = B, n1 = n1[j], n2 = n2[j], seed = i + 2*B)
    
    
    # P-values
    res_Normal_app2[i,j] = (sum(emp_density_Norm_ap >= 0)+1)/(B + 1)
    res_PJ2[i,j] = (sum(emp_density_PJ >= 0)+1)/(B + 1)
    res_NPJ2[i,j] = (sum(emp_density_NPJ >= 0)+1)/(B + 1)
    
    
    x = rep(0, n1[j])
    y = rep(0, n2[j])
    x[runif(n1[j]) < theta1] = 1
    y[runif(n2[j]) < theta2] = 1
    res_np2[i,j] = prop.test(x = c(sum(x), sum(y)), n = c(n1[j], n2[j]), alternative = "less")$p.value
    
  }
  
}

# apply(res2 < 0.05, 2, mean, na.rm = T)
# apply(res_np2 < 0.05, 2, mean, na.rm = T)

# Compute the power
power_res_Normal2 <- apply(res_Normal_app2 < 0.05, 2, mean, na.rm = TRUE)
power_res_PJ2 <- apply(res_PJ2 < 0.05, 2, mean, na.rm = TRUE)
power_res_NPJ2 <- apply(res_NPJ2 < 0.05, 2, mean, na.rm = TRUE)
power_res_np2 <- apply(res_np2 < 0.05, 2, mean, na.rm = TRUE)


# # Define the n1 vector for x-axis
# n1 <- 16*2^(0:5)

# Plot the graph
plot(n1, power_res_Normal2, type = "b", pch = 19, col = "blue", xlab = "sample size", ylab = "Power", ylim = c(0,1), main = "Power of Private and Non-priavte Test (H=10^3, B=10^3, eps=1)")
points(n1, power_res_np2, type = "b", pch = 19, col = "red")
points(n1, power_res_PJ2, type = "b", pch = 19, col = "green")
points(n1, power_res_NPJ2, type = "b", pch = 19, col = "purple")
legend("bottomright", legend = c("Private Normal Approx", "Non-private", "Private Fudicial", "Non-private Fudicial"), col = c("blue", "red", "green", "purple"), pch = 19, cex = 0.8)


# # Record the end time
# end_time <- Sys.time()
# 
# # Calculate the execution time
# execution_time <- end_time - start_time
# print(paste("Execution time:", execution_time))


# ### To Run on HPC
# # Plot file name
# plot_filename <- "power_plot.png"
# 
# # Open PNG device for plotting
# png(plot_filename)
# 
# # Plot the graph
# plot(n1, power_res_Normal2, type = "b", pch = 19, col = "blue", xlab = "sample size", ylab = "Power", ylim = c(0,1), main = "Power of Private and Non-priavte Test (H=10^3, B=10^3, eps=1)")
# points(n1, power_res_np2, type = "b", pch = 19, col = "red")
# points(n1, power_res_PJ2, type = "b", pch = 19, col = "green")
# points(n1, power_res_NPJ2, type = "b", pch = 19, col = "purple")
# legend("bottomright", legend = c("Private Normal Approx", "Non-private", "Private Fudicial", "Non-private Fudicial"), col = c("blue", "red", "green", "purple"), pch = 19, cex = 0.8)
# 
# # Close the PNG device
# dev.off()
