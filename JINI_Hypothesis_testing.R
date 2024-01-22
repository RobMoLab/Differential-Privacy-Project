# Load the required libraries
library(ggplot2)
library(tidyr)
library(LaplacesDemon)

######### Jordan's functions
FU = function(u){# domain is -1/2 to 1/2
  return(u+1/2)
}

ptulap = function(t,b){
  cdf = function(t,b){
    ifelse(t<=0,b^(-round(t))/(1+b)*(b+FU(t-round(t))*(1-b)),1-b^(round(t))/(1+b)*(b+FU(round(t)-t)*(1-b)))
  }
  return(sapply(t,cdf,b=b))
}

######### Indirect functions

## Private Estimate
get_pi = function(theta, u, w, eps, n){
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}

## Indirect Estimate
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

## Distribution of Indirect Estimator
jini = function(pi0, B, eps, n, seed){
  set.seed(seed)
  res = rep(NA, B)
  counter = 0
  
  while (counter < B){
    inter = get_theta_tilde(pi0 = pi0, Zi = rnorm(1), Wj = runif(1, -0.5, 0.5), eps = eps, n = n)
    
    if (inter < 100){
      counter = counter + 1
      res[counter] = inter
    }
  }
  res
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


# # NON-private Fudicial Application Setting
# JINI_algorithm_Fudicial_Nonprivate<-function(pi0, B, n, seed){
#   res = rep(NA, B)
#   for (i in 1:B) {
#     res[i]<-rbeta(1, n*pi0+0.5, n*(1-pi0)+0.5, ncp = 0)
#     
#   }
#   res
# }












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
      JINI_solution[i] <- 0#in(simulated_sample)
    } else {
      JINI_solution[i] <- inter(simulated_sample, pi_j_star)
    }
  }
  JINI_solution
}



#################################################################
 # Power Simulation for the test; H0: theta <= 0.9 vs Ha: theta > 0.9

H = 10^3
B = 10^3 
n = 16*2^(0:5)
theta0 = 0.95
eps = 1
res_jini =res_jini_normal_approx=res_jini_NP=res_jini_inter=res_nonpriv = res_tulap = res_norm = matrix(NA, H, length(n))


for (j in 1:length(n)) {
  
  for (i in 1:H) {
    set.seed(i)
    
    pi0 = get_pi(theta = theta0, u = runif(n[j]), w = runif(1, -0.5, 0.5), eps = eps, n = n[j])
    
    # # Nonprivate JINI P-Value
    # emp_density_Np = JINI_algorithm_Fudicial_Nonprivate(pi0 = theta0, B = B, n = n[j], seed = i + 2*B)
    # res_jini_NP[i,j] = sum(emp_density_Np <= 0.9)/(B+1)
    
    # JINI P-Value
    emp_density = JINI_algorithm_Fudicial_Appl(pi0 = pi0, B = B, eps = eps, n = n[j], seed = i + 2*B)
    res_jini[i,j] = sum(emp_density <= 0.9)/(B+1)
    
    
    # JINI with normal approximation of binomial
    emp_density_normal_approx = jini(pi0 = pi0, B = B, eps = eps, n = n[j], seed = i + 2*B)
    res_jini_normal_approx[i,j] = sum(emp_density_normal_approx <= 0.9)/(B + 1)
    
    
    # JINI with interpolation
    emp_density_inter = JINI_Inter(pi0 = pi0, B = B, eps = eps, n = n[j], seed = i + 2*B)
    res_jini_inter[i,j] = sum(emp_density_inter <= 0.9)/(B + 1)
    
    
    ###### Jordan's Estimators
    
    X = rbinom(n = 1, size = n[j], prob = theta0)
    values = seq(0, n[j])
    pdf = dbinom(values, size = n[j], prob = 0.9)
    
    # Non-Private P-Value
    Z = X + runif(n = 1, min = -1/2, max = 1/2)
    cdf = punif(values-Z, min = -1/2, max = 1/2)
    res_nonpriv[i,j] = t(cdf)%*%pdf
    
    # Tulap P-Value
    U = runif(n = 1, min = -1/2, max = 1/2)
    G1 = rgeom(n = 1, prob = 1 - exp(-eps))
    G2 = rgeom(n = 1, prob = 1 - exp(-eps))
    Tulap = X + U + G1 - G2
    cdf = ptulap(values - Tulap, exp(-eps))
    res_tulap[i,j] = t(cdf)%*%pdf
    
    # Normal-Normal P-Value
    E1 = rexp(n = 1, rate = eps)
    E2 = rexp(n = 1, rate = eps)
    Laplace = X + E1 - E2
    res_norm[i,j] = 1 - pnorm(Laplace, m = 0.9*n[j], s = sqrt(n[j]*0.9*(1 - 0.9)+2/eps^2))
    
    
  }
  
}


# Calculate the power for each method
power_jini <- apply(res_jini < 0.05, 2, mean)
power_jini_normal_approx <- apply(res_jini_normal_approx < 0.05, 2, mean)
#power_jini_NP <- apply(res_jini_NP < 0.05, 2, mean)
power_jini_inter <- apply(res_jini_inter < 0.05, 2, mean)
power_nonpriv <- apply(res_nonpriv < 0.05, 2, mean)
power_tulap <- apply(res_tulap < 0.05, 2, mean)
power_norm <- apply(res_norm < 0.05, 2, mean)

# Combine the power calculations with the corresponding sample sizes
power_data <- data.frame(
  SampleSize = rep(n, 6),
  Power = c(power_jini, power_jini_normal_approx,power_jini_inter, power_nonpriv, power_tulap, power_norm),
  Method = factor(rep(c("JINI Fudicial", "JINI with Normal approximation to Binomial", "JINI Interpolation",  "Non-Private", "Tulap", "Normal-Normal"), each = length(n)))
)



ggplot(power_data, aes(x = SampleSize, y = Power, color = Method)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(trans = 'log2', breaks = n) +  # Assuming n is in powers of 2
  labs(title = "Power (B=H=10^3, eps=1, theta0 = 0.95)",
       x = "Sample Size (n)",
       y = "Power") +
  theme_minimal() +
  theme(legend.position = "bottom")







##########################################################################################################################
# Level Simulation for the test; H0: theta <= 0.9 vs Ha: theta > 0.9

H = 10^3
B = 10^3
n = 30 
theta0 = seq(0.1, 0.9, by = 0.1)
eps = 1

res_jini =res_jini_normal_approx=res_jini_NP=res_jini_inter=res_nonpriv = res_tulap = res_norm = matrix(NA, H, length(theta0))


for (j in 1:length(theta0)){
  
  
  for (i in 1:H){
    
    set.seed(i)
    pi0 = get_pi(theta = theta0[j], u = runif(n), w = runif(1, -0.5, 0.5), eps = eps, n = n)
    
    # JINI P-Value
    emp_density = JINI_algorithm_Fudicial_Appl(pi0 = pi0, B = B, eps = eps, n = n, seed = i + 2*B)
    res_jini[i,j] = sum(emp_density <= theta0[j])/(B+1)
    
    # # Nonprivate JINI P-value
    # emp_density_JINI_NP = JINI_algorithm_Fudicial_Nonprivate(pi0 = pi0, B = B, n = n, seed = i + 2*B)
    # res_jini_NP[i,j] = sum(emp_density_JINI_NP <= theta0[j])/(B+1)
    
    
    # JINI with normal approximation to binomial P-Value  
    emp_density_normal_approx = jini(pi0 = pi0, B = B, eps = eps, n = n, seed = i + 2*B)
    res_jini_normal_approx[i,j] = sum(emp_density_normal_approx <= theta0[j])/(B + 1)
    
  
    # JINI with interpolation
    emp_density_inter = JINI_Inter(pi0 = pi0, B = B, eps = eps, n = n, seed = i + 2*B)
    res_jini_inter[i,j] = sum(emp_density_inter <= theta0[j])/(B + 1)
    
    
    
    ###### Jordan's Estimators
    
    X = rbinom(n = 1, size = n, prob = theta0[j])
    values = seq(0, n)
    pdf = dbinom(values, size = n, prob = theta0[j])
    
    # Non-Private P-Value
    Z = X + runif(n = 1, min = -1/2, max = 1/2)
    cdf = punif(values-Z, min = -1/2, max = 1/2)
    res_nonpriv[i,j] = t(cdf)%*%pdf
    
    # Tulap P-Value
    U = runif(n = 1, min = -1/2, max = 1/2)
    G1 = rgeom(n = 1, prob = 1 - exp(-eps))
    G2 = rgeom(n = 1, prob = 1 - exp(-eps))
    Tulap = X + U + G1 - G2
    cdf = ptulap(values - Tulap, exp(-eps))
    res_tulap[i,j] = t(cdf)%*%pdf
    
    # Normal-Normal P-Value
    E1 = rexp(n = 1, rate = eps)
    E2 = rexp(n = 1, rate = eps)
    Laplace = X + E1 - E2
    res_norm[i,j] = 1 - pnorm(Laplace, m = theta0[j]*n, s = sqrt(n*theta0[j]*(1 - theta0[j])+2/eps^2))
    
    print(i/H)
    
  }
  
  print(j/length(theta0))
  
}


# Calculate the level for each method
type1_error_jini <- apply(res_jini < 0.05, 2, mean)
#type1_error_jini_NP <- apply(res_jini_NP < 0.05, 2, mean)
type1_error_jini_normal_approx <- apply(res_jini_normal_approx < 0.05, 2, mean)
type1_error_jini_inter <- apply(res_jini_inter < 0.05, 2, mean)
type1_error_nonpriv <- apply(res_nonpriv < 0.05, 2, mean)
type1_error_tulap <- apply(res_tulap < 0.05, 2, mean)
type1_error_norm <- apply(res_norm < 0.05, 2, mean)


# data frame of level
level_data <- data.frame(
  Theta0 = rep(theta0, each = 6),
  Type1Error = c(type1_error_jini, type1_error_jini_normal_approx, type1_error_jini_inter,
                 type1_error_nonpriv, type1_error_tulap, type1_error_norm),
  Method = rep(c("JINI Fudicial", "JINI Normal Approx to Normal", "JINI Interpolation", "Non-Private", "Tulap", "Normal-Normal"), 
               times = length(theta0))
)





ggplot(level_data, aes(x = Theta0, y = Type1Error, color = Method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "orange") +
  scale_y_continuous(limits = c(0.01, 0.07)) +
  labs(title = "Type I Error Rate  (H=B=10^3, eps=1, n=30)",
       x = "Theta0",
       y = "Type I Error Rate") +
  theme_minimal() +
  theme(legend.position = "bottom")




