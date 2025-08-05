library(binomialDP)
library(dplyr)
library(progress)



### Helping function ####
FU <- function(u) {
  return(u + 1/2)
}

#' Tulap Cumulative Distribution Function
#' 
#' Computes the CDF of the Tulap (Truncated-Uniform-Laplace) distribution.
#' 
#' @param t Evaluation points
#' @param b Parameter controlling noise level (b = exp(-ε))
#' @return Vector of CDF values
ptulap = function(t,b){
  cdf = function(t,b){
    ifelse(t<=0,b^(-round(t))/(1+b)*(b+FU(t-round(t))*(1-b)),1-b^(round(t))/(1+b)*(b+FU(round(t)-t)*(1-b)))
  }
  return(sapply(t,cdf,b=b))
}




get_pi <- function(theta, eps, n) {
  # Input validation
  if (theta <= 0 || theta >= 1) stop("theta must be between 0 and 1")
  if (eps <= 0) stop("eps must be positive")
  if (n <= 0) stop("n must be positive")
  
  w <- runif(1, -0.5, 0.5)
  u <- runif(n)
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}



#'FIMA
#'
#' Generates fiducial samples for differentially private inference using
#' a fiducial approach with beta distributions.
#'
#' @param pi0 Private proportion estimate
#' @param eps Privacy parameter (ε)
#' @param n Sample size
#' @param H Number of fiducial samples (default = 1)
#' @param delta Sensitivity parameter (default = 1)
#' @param seed Random seed (default = 123)
#'
#' @return Vector of H fiducial samples
#' @examples
#' fima(pi0 = 0.5, eps = 1, n = 100, H = 1000)
fima <- function(pi0, eps, n, H , delta = 1, seed = 123) {
  
  set.seed(seed)
  
  # Generate noise variables
  Wj <- runif(H, -0.5, 0.5)
  
  # Create perturbed proportions
  pi_star <- pi0 + delta/(eps*n)*sign(Wj)*log(1 - 2*abs(Wj))
  
  # Initialize fiducial samples
  theta_ifb <- rep(NA, H)
  
  # Generate beta samples for valid proportions
  index <- pi_star < 1 & pi_star > 0
  theta_ifb[index] <- rbeta(
    sum(index), 
    n*pi_star[index] + 0.5,  # Add 0.5 for Haldane-Anscombe correction
    n*(1 - pi_star[index]) + 0.5
  )
  
  # Handle boundary cases
  theta_ifb[pi_star >= 1] <- 1 - .Machine$double.eps
  theta_ifb[pi_star <= 0] <- .Machine$double.eps
  
  return(theta_ifb)
}




simulate_all_methods <- function(theta0, n, B, eps, H, alpha = 0.05) {
  # Initialize results storage
  results <- expand.grid(
    theta = theta0,
    rep = 1:B,
    method = c("FIMA", "NonPrivate", "Tulap", "Exact"),
    type1_error = NA,
    time_ms = NA,
    stringsAsFactors = FALSE
  )
  
  pb <- progress_bar$new(total = length(theta0) * B, format = "[:bar] :percent")
  
  for (j in seq_along(theta0)) {
    for (i in 1:B) {
      set.seed(i)
      
      # Generate data
      x <- rbinom(1, size = n, prob = theta0[j])
      values <- seq(0, n)
      pdf <- dbinom(values, size = n, prob = theta0[j])
      
      # --- FIMA Method ---
      start_time <- Sys.time()
      pi0 <- get_pi(theta = theta0[j], eps = eps, n = n)
      fima_samples <- fima(pi0 = pi0, H = H, eps = eps, n = n, seed = i)
      fima_pval <- mean(fima_samples > theta0[j])
      fima_time <- as.numeric(Sys.time() - start_time) * 1000
      
      # --- Non-Private Method ---
      start_time <- Sys.time()
      Z <- x + runif(1, -1/2, 1/2)
      np_pval <- sum(punif(values - Z, min = -1/2, max = 1/2) * pdf)
      np_time <- as.numeric(Sys.time() - start_time) * 1000
      
      # --- Tulap Method ---
      start_time <- Sys.time()
      U <- runif(1, -1/2, 1/2)
      G1 <- rgeom(1, prob = 1 - exp(-eps))
      G2 <- rgeom(1, prob = 1 - exp(-eps))
      Tulap <- x + U + G1 - G2
      cdf <- ptulap(values - Tulap, exp(-eps))
      tulap_pval <- sum(cdf * pdf)  
      tulap_time <- as.numeric(Sys.time() - start_time) * 1000
      
      # --- Exact Method ---
      start_time <- Sys.time()
      exact_pval <- binom.test(x, n, p = theta0[j])$p.value
      exact_time <- as.numeric(Sys.time() - start_time) * 1000
      
      # Store results
      idx <- which(results$theta == theta0[j] & results$rep == i)
      results[idx, ] <- data.frame(
        theta = theta0[j],
        rep = i,
        method = c("FIMA", "NonPrivate", "Tulap", "Exact"),
        type1_error = c(
          fima_pval < alpha,
          np_pval < alpha,
          tulap_pval < alpha,
          exact_pval < alpha
        ),
        time_ms = c(fima_time, np_time, tulap_time, exact_time)
      )
      
      pb$tick()
    }
  }
  return(results)
}

# Run simulation 
theta0 <- seq(0.1, 0.9, by = 0.1)
n <- 30
B <- 10^4  
eps <- 1
H <- 10^3

results <- simulate_all_methods(theta0, n, B, eps, H)

# Analyze results


# Proper summary statistics with theta breakdown
full_stats <- results %>%
  group_by(theta, method) %>%
  summarise(
    type1_rate = mean(type1_error),
    type1_se = sd(type1_error)/sqrt(n()),
    avg_time = mean(time_ms),
    time_se = sd(time_ms)/sqrt(n()),
    iter = n(),
    .groups = "drop"
  )


NonPrivate_stats <- full_stats  %>% filter(method == "NonPrivate")
exact_stats <- full_stats  %>% filter(method == "Exact")
tulap_stats <- full_stats  %>% filter(method == "Tulap")
fima_stats <- full_stats  %>% filter(method == "FIMA")
