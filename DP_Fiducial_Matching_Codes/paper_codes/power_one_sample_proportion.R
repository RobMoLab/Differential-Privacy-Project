library(dplyr)
library(progress)





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


#' Compute Differentially Private Estimate of Proportion

get_pi <- function(theta, u, w, eps, n) {
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}

#' FIMA algorithm
#'
#' Generates samples from the fiducial distribution for differentially private
#' binomial proportion estimation.
#'
#' @param pi0 Initial private estimate
#' @param eps Privacy parameter (epsilon > 0)
#' @param n Sample size (positive integer)
#' @param H Number of Monte Carlo samples (positive integer, default=1)
#' @param seed Random seed (default=123)
#' @return Vector of H fiducial samples

fima <- function(pi0, eps, n, H = 1, delta = 1, seed = 123) {
  set.seed(seed)
  
  Wj <- runif(H, -0.5, 0.5)
  pi_star <- pi0 + delta/(eps*n)*sign(Wj)*log(1 - 2*abs(Wj))
  
  theta_ifb <- rep(NA, H)
  
  index <- pi_star < 1 & pi_star > 0
  theta_ifb[index] <- rbeta(length(index), 
                            n*pi_star[index] + 0.5, 
                            n*(1 - pi_star[index]) + 0.5)
  theta_ifb[pi_star >= 1] <- 1 - .Machine$double.eps
  theta_ifb[pi_star <= 0] <- .Machine$double.eps
  
  return(theta_ifb)
}




run_combined_simulation <- function(theta0, theta_values, n, B, eps = 1, H = 1000, alpha = 0.05) {
  # Initialize results storage
  results <- expand.grid(
    theta0 = theta0,
    theta = theta_values,
    rep = 1:B,
    method = c("Exact", "Tulap", "NonPrivate", "FIMA"),
    rejected = NA,
    time_ms = NA,
    stringsAsFactors = FALSE
  )
  
  pb <- progress_bar$new(total = length(theta_values) * B, format = "[:bar] :percent")
  
  for (j in seq_along(theta_values)) {
    for (i in 1:B) {
      set.seed(i)
      
      # Generate data per iteration
      x <- rbinom(1, size = n, prob = theta_values[j])
      values <- seq(0, n)
      pdf_H0 <- dbinom(values, size = n, prob = theta0)
      
      # --- Exact Binomial Test ---
      start_time <- Sys.time()
      exact_pval <- binom.test(x, n, p = theta0, alternative = "greater")$p.value
      exact_time <- as.numeric(Sys.time() - start_time) * 1000
      
      # --- Tulap Method ---
      start_time <- Sys.time()
      U <- runif(1, -0.5, 0.5)
      G1 <- rgeom(1, prob = 1 - exp(-eps))
      G2 <- rgeom(1, prob = 1 - exp(-eps))
      noisy_x <- x + U + G1 - G2
      cdf <- ptulap(values - noisy_x, exp(-eps))
      tulap_pval <- sum(cdf * pdf_H0)
      tulap_time <- as.numeric(Sys.time() - start_time) * 1000
      
      # --- Non-Private Method ---
      start_time <- Sys.time()
      Z <- x + runif(1, -0.5, 0.5)
      np_pval <- sum(punif(values - Z, -0.5, 0.5) * pdf_H0)
      np_time <- as.numeric(Sys.time() - start_time) * 1000
      
      # --- FIMA Method ---
      start_time <- Sys.time()
      pi0 <- get_pi(theta = theta_values[j], u = runif(n), w = runif(1, -0.5, 0.5), eps = eps, n = n)
      fima_samples <- fima(pi0 = pi0, H = H, eps = eps, n = n, seed = i + 2*H)
      fima_pval <- mean(fima_samples <= theta0)
      fima_time <- as.numeric(Sys.time() - start_time) * 1000
      
      # Store results
      idx <- which(results$theta == theta_values[j] & results$rep == i)
      results[idx, ] <- data.frame(
        theta0 = theta0,
        theta = theta_values[j],
        rep = i,
        method = c("Exact", "Tulap", "NonPrivate", "FIMA"),
        rejected = c(
          exact_pval < alpha,
          tulap_pval < alpha,
          np_pval < alpha,
          fima_pval < alpha
        ),
        time_ms = c(exact_time, tulap_time, np_time, fima_time)
      )
      
      pb$tick()
    }
  }
  
  return(results)
}

# Analysis Functions
analyze_results <- function(sim_results) {
  # Overall summary
  overall <- sim_results %>%
    group_by(method) %>%
    summarise(
      power = mean(rejected),
      power_se = sd(rejected)/sqrt(n()),
      avg_time_ms = mean(time_ms),
      time_se = sd(time_ms)/sqrt(n()),
      .groups = "drop"
    )
  
  # Breakdown by theta
  by_theta <- sim_results %>%
    group_by(theta, method) %>%
    summarise(
      power = mean(rejected),
      power_se = sd(rejected)/sqrt(n()),
      avg_time_ms = mean(time_ms),
      .groups = "drop"
    )
  
  return(list(overall = overall, by_theta = by_theta))
}



theta0 <- 0.2
theta_values <- seq(0.2, 0.9, by = 0.1)
n <- 30
B <- 10^4  
eps <- 1
H <- 10^3

results <- run_combined_simulation(theta0, theta_values, n, B, eps, H)





# Full analysis 
full_results <- results %>%
  group_by(theta, method) %>%
  summarise(
    power = mean(rejected),
    power_se = sd(rejected)/sqrt(n()),
    mean_time_ms = mean(time_ms),
    time_se = sd(time_ms)/sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Create filtered datasets for each method
exact_results <- full_results %>% filter(method == "Exact")
tulap_results <- full_results %>% filter(method == "Tulap")
np_results <- full_results %>% filter(method == "NonPrivate")
fima_results <- full_results %>% filter(method == "FIMA")
