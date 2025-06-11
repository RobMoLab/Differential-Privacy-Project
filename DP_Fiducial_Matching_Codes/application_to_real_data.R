
#####Helping Functions######
# Function to compute chi2 stat. based on counts in contingency table
compute_chi2 <- function(contingency_tab){
  n <- sum(contingency_tab)
  marginal_row <- apply(contingency_tab,1,sum)/n 
  maringal_col <- apply(contingency_tab,2,sum)/n
  expected <- outer(marginal_row,maringal_col,FUN="*")*n
  return(sum((contingency_tab-expected)^2/expected))
}
# Laplace noise
laplace_noise <- function(sensitivity, epsilon){
  W <- runif(1, -0.5, 0.5)
  noise <- -sensitivity/epsilon*sign(W)*log(1 - 2*abs(W))
  return(noise)
}
# Laplace mechanism for differential privacy
add_laplace_noise <- function(counts, epsilon) {
  sensitivity <- 2
  noisy_counts <- sapply(counts, function(count) count+laplace_noise(sensitivity, epsilon=epsilon))
  noisy_counts <- pmax(noisy_counts, 0)
  return(haldane_anscombe_correction(noisy_counts))  # Ensure non-negative counts
}
# Haldane anscombe correction for when cells are equal to 0
haldane_anscombe_correction <- function(counts){
  counts[counts==0] <- counts[counts==0]+0.5
  return(counts)
}


# get private proportion
get_pi <- function(theta_hat, eps, n) {
  # Generate a value for w until it satisfies the condition
  repeat {
    w <- runif(1, -1, 1)
    if (1 - 2 * abs(w) > 0) {
      break
    }
  }
  
  # obtain laplace noise to yield epsilon-DP
  noise <- -1 / (eps * n) * sign(w) * log(1 - 2 * abs(w))
  theta_hat + noise
}


#pi0<-get_pi(0.6, u=runif(30), w=runif(1,-1,1), eps = 1, n=30)


# Private
JINI_algorithm_Fiducial_Appl <-function(pi0, H, eps, n, seed){
  set.seed(seed)
  res = rep(NA, H)
  counter = 0
  
  while (counter < H){
    #set.seed(seed *(1+counter))
    Wj = runif(1, -0.5, 0.5)
    pi_j_star = pi0 +1/(eps*n)*sign(Wj)*log(1 - 2*abs(Wj))
    
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

jini2samples_fiducial_Priv = function(pi1, pi2, H, eps, n1, n2){
  #set.seed(seed)
  res = rep(NA, H)
  
  distri_grp1 = JINI_algorithm_Fiducial_Appl(pi0 = pi1, H=H, eps = eps/2, n = n1, seed=sample(1:H,1))
  distri_grp2 = JINI_algorithm_Fiducial_Appl(pi0 = pi2, H=H, eps = eps/2, n = n2, seed=sample(H+1:2*H,1))
  
  
  res = distri_grp1 - distri_grp2
  return(res)
}



#####     HIV diagnoses in the US and 6 territories and freely associated states for the most-affected subpopulations, 2022  #####
######                https://www.cdc.gov/hiv/data-research/facts-stats/index.html                                             #####

# H_0:p_1=p_2 vs H_1: p_1>p_2
# Parameters/data
eps =0.1
total_infected=31800
Hispanic_Latino = 9374
Black_Afri_Ameri=8831

p1_hat = Hispanic_Latino/total_infected
p2_hat = Black_Afri_Ameri/total_infected
H=10000
n=total_infected

# Private proportions
pi1 =get_pi( p1_hat, eps/2, n)
pi2 =get_pi(p2_hat, eps/2, n)

# Fima p-values
p_value<-mean(jini2samples_fiducial_Priv(pi1, pi2, H,eps,n,n) <=0 )
p_value



### Nonprivate tests  ########

# Nonprivate z-test

# Calculate pooled proportion
p_pooled <- (p1_hat * n + p2_hat * n) / (n + n)

# Calculate standard error
SE <- sqrt(p_pooled * (1 - p_pooled) * (1/n + 1/n))

# Calculate z-score (test statistic)
z_score <- (p1_hat - p2_hat) / SE

# Calculate p-value
p_value <- 1-pnorm(z_score)

# Print results
#cat("Z-score:", z_score, "\n")
cat("P-value:", p_value, "\n")


  ### Exact Binomial test
# Calculate the number of successes and failures
x1 <- round(p1_hat * n)  # Number of successes in group 1
x2 <- round(p2_hat * n)  # Number of successes in group 2
failure1 <- n - x1   # Number of failures in group 1
failure2 <- n - x2   # Number of failures in group 2

# Create contingency table
table <- matrix(c(x1, x2, failure1, failure2), nrow = 2, byrow = TRUE)

# Perform Fisher's exact test
fisher_test <- fisher.test(table, alternative = "greater")

non_private_p_value<- fisher_test$p.value

# Print results
cat("P-value:", non_private_p_value, "\n")


### Evealuating consistency ####

multiple_trials_two_sample <- function(p1_hat, p2_hat, H,eps,n1,n2, num_run){
  p_value_run <- c()
  seeds<-sample(1:num_run,num_run)
  for (i in 1:num_run) {
    set.seed(i)
    pi1 =get_pi( p1_hat, eps/2, n)
    pi2 =get_pi(p2_hat, eps/2, n)
    p_value_run[i]<-mean(jini2samples_fiducial_Priv(pi1, pi2, H,eps,n1,n2) <=0)

  }
  #ave_p_value<-mean(p_value_run<=0.002)
  min_p_value = min(p_value_run)
  max_p_value = max(p_value_run)
  cat("min p_value:", min_p_value, "\nmax p_value:", max_p_value)

}

multiple_trials_two_sample(pi1, pi2, H,eps,n,n,10000)



# Effect of epsilon
two_sample_diff_epsilon<-function(pi1, pi2, H,eps_values,n1,n2){
  k=length(eps_values)
  p_value_eps <- c()
  set.seed(123)
  for (i in 1:k) {
    eps= eps_values[i]
    pi1 =get_pi( p1_hat, eps/2, n)
    pi2 =get_pi(p2_hat, eps/2, n)
    p_value_eps[i]<-mean(jini2samples_fiducial_Priv(pi1, pi2, H,eps,n1,n2) <=0)
  }
  return(p_value_eps)
}

eps_values=c(0.001, 0.01, 0.1,1, 0.5, 3, 5, 10 )
two_sample_diff_epsilon(pi1, pi2, H,eps_values,n,n)





                               ################### Alabama HIV data  ##########################

# H_0:p=0.7  Vs H_1: p<0.7

specified_value = 0.7
observed_count =232
total = 374
epsilon=1
H=10000
seed=123
alpha <- 0.05   
oberseved_prop =observed_count/total
set.seed(seed)


# private test with Fima
private_observed_prop = get_pi(oberseved_prop, epsilon,total )


distri = JINI_algorithm_Fiducial_Appl(private_observed_prop, H, epsilon, total, seed)

private_p_value = mean(distri>=0.7)
private_p_value


        # Nonprivate t-test
# Calculate test statistic
se <- sqrt((specified_value * (1 - specified_value)) / total)
t_statistic <- (oberseved_prop - specified_value) / se

# Calculate p-value for left-sided test
p_value <- pt(t_statistic, df = total - 1, lower.tail = TRUE)

# Output results
cat("P-value:", p_value, "\n")



# Private CI
quantile(distri, probs = c(alpha/2, 1-alpha/2))


# Nonprivate CI
proportion_conf_interval <- function(observed_proportion, n, alpha = 0.05) {
  # Standard error for the proportion
  standard_error <- sqrt(observed_proportion * (1 - observed_proportion) / n)
  
  # Z-score for the given alpha level (two-sided)
  z_score <- qnorm(1 - alpha / 2)
  
  # Confidence interval bounds
  lower_bound <- observed_proportion - z_score * standard_error
  upper_bound <- observed_proportion + z_score * standard_error
  
  return(c(lower_bound, upper_bound))
}

proportion_conf_interval(oberseved_prop, total, alpha)







#### Is race associated with poverty status?

#RACE/ORIGIN                                         Below poverty level      Not below poverty level      Total

# White alone                                         21,525,577                                            213,295,033
# Black or African American alone                     8,519,391                                             39,695,427
# American Indian and Alaska Native alone             608,547                                               2,692,978
# Asian alone                                         1,897,150                                             18,754,209
# Native Hawaiian and Other Pacific Islander alone    103,050                                               607,291
# Some other race alone                               3,652,060                                             19,671,062
# Two or more races                                   4,215,809                                             28,559,448
# Hispanic or Latino origin (of any race)             10,447,540                                            60,614,309
# White alone, not Hispanic or Latino                 17,620,793                                            190,513,343


non_private_table_US<-matrix(c(21525577, 8519391,  608547,  1897150, 103050, 3652060, 4215809, 10447540, 17620793,
                               213295033,  39695427, 2692978, 18754209, 607291, 19671062, 28559448,   60614309, 190513343), nrow = 9)

below_poverty_level =non_private_table_US[,2]-non_private_table_US[,1]

non_private_table_US[,2]=below_poverty_level

#data <- read_excel("/home/ogonna/Differential_Privacy_Project/Privacy_Project/Applications_to_real_data/ACSST5Y2022.S1701-2024-11-01T171755.xlsx", sheet = "Data")

#data[colnames(data)=="United States"]


# private chi-square test with Fima
Private_chisquare_test <- function(counts, epsilon, H) {
  
  chi2_hb <- rep(NA, H)
  noisy_contingency_table <- add_laplace_noise(as.vector(counts), epsilon)
  noisy_contingency_table<-matrix( noisy_contingency_table, nrow = dim(counts)[1])
  chi2_0b <- compute_chi2( noisy_contingency_table)
  
  for (h in 1:H) {
    denoised_counts <- add_laplace_noise(noisy_contingency_table, epsilon)
    n_star <- sum(denoised_counts)
    denoised_counts<-matrix(denoised_counts, nrow =dim(counts)[1] )
    marginal_row <- apply(denoised_counts, 1, sum) / n_star
    marginal_col <- apply(denoised_counts, 2, sum) / n_star  
    
    marginal_row_tilde <- sapply(as.vector(marginal_row), function(p) rbeta(1, n_star*p + 0.5, n_star*(1 - p) + 0.5, ncp = 0))
    marginal_col_tilde <- sapply(as.vector(marginal_col), function(p) rbeta(1, n_star*p + 0.5, n_star*(1 - p) + 0.5, ncp = 0))
    
    joint_probs0_tilde <- outer(marginal_row_tilde, marginal_col_tilde)
    counts_tilde_star <- rmultinom(1, n_star, prob = joint_probs0_tilde)
    counts_tilde_star <- haldane_anscombe_correction(counts_tilde_star)
    counts_tilde <- matrix(add_laplace_noise(counts_tilde_star, epsilon), nrow = dim(counts)[1], byrow = F)
    chi2_hb[h] <- compute_chi2(counts_tilde)
  }
  p_value = mean(chi2_hb >= chi2_0b)
  
  #cat("p-value:", p_value)
  return(p_value)
  
}


Private_chisquare_test(non_private_table_US,0.00001,1000)


# Nonprivate chi-square test
chisq.test(non_private_table_US)$p.value



# Effect of epsilon
Private_chisquare_test_eps <- function(counts, eps_values, H) {
  p_values_eps <- numeric(length(eps_values))
  for(i in seq_along(eps_values)) {
    epsilon <- eps_values[i]
    p_values_eps[i] <- Private_chisquare_test(counts, epsilon, H)
  }
  return(p_values_eps)
}

eps_values=c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5,1,3,5,10)
H=10^3
Private_chisquare_test_eps(non_private_table_US,eps_values,H)




# Evaluating the consistency of the test
Private_chisquare_test_trial <- function(counts, epsilon, H, n_run) {
  p_value_run=c()
  chi2_hb <- rep(NA, H)
  seeds = sample(1:n_run, n_run)
  for (i in 1:n_run) {
    set.seed(seeds[i])
    noisy_contingency_table <- add_laplace_noise(as.vector(counts), epsilon)
    noisy_contingency_table<-matrix( noisy_contingency_table, nrow = dim(counts)[1])
    chi2_0b <- compute_chi2( noisy_contingency_table)
    
    for (h in 1:H) {
      denoised_counts <- add_laplace_noise(noisy_contingency_table, epsilon)
      n_star <- sum(denoised_counts)
      denoised_counts<-matrix(denoised_counts, nrow =dim(counts)[1] )
      marginal_row <- apply(denoised_counts, 1, sum) / n_star
      marginal_col <- apply(denoised_counts, 2, sum) / n_star  
      
      marginal_row_tilde <- sapply(as.vector(marginal_row), function(p) rbeta(1, n_star*p + 0.5, n_star*(1 - p) + 0.5, ncp = 0))
      marginal_col_tilde <- sapply(as.vector(marginal_col), function(p) rbeta(1, n_star*p + 0.5, n_star*(1 - p) + 0.5, ncp = 0))
      
      joint_probs0_tilde <- outer(marginal_row_tilde, marginal_col_tilde)
      counts_tilde_star <- rmultinom(1, n_star, prob = joint_probs0_tilde)
      counts_tilde_star <- haldane_anscombe_correction(counts_tilde_star)
      counts_tilde <- matrix(add_laplace_noise(counts_tilde_star, epsilon), nrow = dim(counts)[1], byrow = F)
      chi2_hb[h] <- compute_chi2(counts_tilde)
    }
    p_value_run[i] = mean(chi2_hb >= chi2_0b)
    
    
  }
  max_pvalue = max(p_value_run)
  min_pvalue = min(p_value_run)
  cat("min p_value:", min_pvalue, "\nmax p_value:", max_pvalue)
  #ave_p_value=mean(p_value_run)
  #print(ave_p_value)
  #return(p_value_run)
}

Private_chisquare_test_trial(non_private_table_US,0.00001,1000,1000)  




                        ###### Alabama poverty level data #######
non_private_table_AL<-matrix(c(3254319, 1273760,20496, 68150, 1999, 92356, 179347,225910, 3175893, 
                               369687, 322862, 3832, 8370, 530,27419,36197,60468,351213), nrow=9, ncol=2)  

Private_chisquare_test(non_private_table_AL, 0.001,1000)

chisq.test(non_private_table_AL)$p.value





                              ### Auburn, AL 36830   Poverty level data #####

non_private_table_zip_code_36830_AL <- matrix(c(6542,1056,13,475,0,41,212,118,6512, 33791, 6516, 51, 4181, 0, 385, 1478, 1286, 33391),
                                              nrow = 9, ncol=2)
non_private_table_zip_code_36830_AL[,2]<-non_private_table_zip_code_36830_AL[,2]-non_private_table_zip_code_36830_AL[,1]
#non_private_table_zip_code_36830_AL

#rowSums(non_private_table_zip_code_36830_AL)

Private_chisquare_test(non_private_table_zip_code_36830_AL, 1,1000)
chisq.test(non_private_table_zip_code_36830_AL)$p.value






