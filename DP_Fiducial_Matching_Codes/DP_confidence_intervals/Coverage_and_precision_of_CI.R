library(binomialDP)
library(progress)
library(dplyr)




# Compute Differentially Private CI Using Tulap Mechanism

tulap_ci <- function(theta0, n, eps, alpha) {
  de <- 0.0                    # Delta parameter for Tulap
  b <- exp(-eps)               # Tulap parameter
  q <- 2 * de * b / (1 - b + 2 * de * b)  # Tulap parameter
  seed <- runif(1, 0, 1)
  N <- rtulap(n = 1, m = 0, b = b, q = q)
  Z <- sdp_fun(seed, N, eps, theta0, n)
  ci <- CITwoSide(alpha = alpha, Z, size = n, b = b, q = q)
  return(c(ci[1], ci[2]))
}

# Helper Function for Binomial SDP Mechanism

sdp_fun <- function(seed, N, ep, theta, n) {
  B <- qbinom(seed, size = n, prob = theta)
  s1 <- B + (1/ep) * N[1]
  return(s1)
}






# Compute Private Estimate of Proportion

get_pi <- function(theta, eps, n) {
  # Input validation
  if (theta <= 0 || theta >= 1) stop("theta must be between 0 and 1")
  if (eps <= 0) stop("eps must be positive")
  if (n <= 0) stop("n must be positive")
  
  w <- runif(1, -0.5, 0.5)
  u <- runif(n)
  mean(u < theta) - 1/(eps*n)*sign(w)*log(1 - 2*abs(w))
}

# FIMA algorithm 

fima <- function(pi0, eps, n, H = 1, delta = 1, seed = 123) {
  # Input validation
  if (eps <= 0) stop("eps must be positive")
  if (n <= 0) stop("n must be positive")
  if (H <= 0) stop("H must be positive")
  
  set.seed(seed)
  Wj <- runif(H, -0.5, 0.5)
  pi_star <- pi0 + delta/(eps*n)*sign(Wj)*log(1 - 2*abs(Wj))
  
  theta_ifb <- rep(NA, H)
  index <- pi_star < 1 & pi_star > 0
  
  tryCatch({
    theta_ifb[index] <- rbeta(
      sum(index), 
      n*pi_star[index] + 0.5, 
      n*(1 - pi_star[index]) + 0.5
    )
  }, error = function(e) {
    stop("Error in beta sampling: ", e$message)
  })
  
  # Handle boundary cases with machine precision
  theta_ifb[pi_star >= 1] <- 1 - .Machine$double.eps
  theta_ifb[pi_star <= 0] <- .Machine$double.eps
  
  return(theta_ifb)
}

# Compute CI for One-sample Proportion Using Fiducial Approach

One_sample_CI_Fiducial <- function(pi0, H, eps, n, alpha, seed) {
  # Input validation
  if (H <= 0) stop("H must be positive")
  if (eps <= 0) stop("eps must be positive")
  if (n <= 0) stop("n must be positive")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  
  tryCatch({
    JINI_solution <- fima(pi0 = pi0, H = H, eps = eps, n = n, seed = seed)
    c(
      quantile(JINI_solution, alpha/2, na.rm = TRUE), 
      quantile(JINI_solution, 1 - (alpha/2), na.rm = TRUE)
    )
  }, error = function(e) {
    warning("Error in CI computation: ", e$message)
    return(c(NA, NA))
  })
}








# Parameters
theta_values <- seq(0.1, 0.98, by = 0.01)  # True proportions
n <- 30                                     # Sample size
B <- 10^4                                  # Replications per theta
alpha <- 0.05                               # Significance level
eps <- 1                                    # Privacy budget (for DP methods)
H <- 10^3                                   # FIMA samples
methods <- c("Wald", "Exact", "Tulap", "FIMA")

# Pre-allocate results (rows = theta_values × B × 4 methods)
results <- expand.grid(
  theta = theta_values,
  rep = 1:B,
  method = methods,
  coverage = NA,
  ci_length = NA,
  time_ms = NA
)

# Progress bar
pb <- progress_bar$new(total = length(theta_values) * B, format = "[:bar] :percent")

# Simulation loop
for (theta in theta_values) {
  for (i in 1:B) {
    set.seed(i)  
    
    # Generate data (x successes in n trials)
    x <- rbinom(1, n, theta)
    
    # (1) Non-private Wald CI (NP in paper)
    start_time <- Sys.time()
    p_hat <- x / n
    se <- sqrt(p_hat * (1 - p_hat) / n)
    margin <- qnorm(1 - alpha / 2) * se
    ci_wald <- c(p_hat - margin, p_hat + margin)
    time_wald <- as.numeric(Sys.time() - start_time) * 1000
    
    # (2) Exact Binomial CI 
    start_time <- Sys.time()
    ci_exact <- binom.test(x, n, conf.level = 1 - alpha)$conf.int
    time_exact <- as.numeric(Sys.time() - start_time) * 1000
    
    # (3) Tulap DP CI
    start_time <- Sys.time()
    ci_tulap <- tulap_ci(theta0 = theta, n = n, eps = eps, alpha = alpha)
    time_tulap <- as.numeric(Sys.time() - start_time) * 1000
    
    # (4) FIMA DP CI
    start_time <- Sys.time()
    pi0 <- get_pi(theta = theta, eps = eps, n = n)
    ci_fima <- One_sample_CI_Fiducial(pi0, H, eps, n, alpha, seed = i)
    time_fima <- as.numeric(Sys.time() - start_time) * 1000
    
    # Store results (find rows for current theta and rep)
    row_idx <- which(results$theta == theta & results$rep == i)
    results[row_idx, ] <- data.frame(
      theta = theta,
      rep = i,
      method = methods,
      coverage = c(
        ci_wald[1] <= theta && theta <= ci_wald[2],
        ci_exact[1] <= theta && theta <= ci_exact[2],
        ci_tulap[1] <= theta && theta <= ci_tulap[2],
        ci_fima[1] <= theta && theta <= ci_fima[2]
      ),
      ci_length = c(
        ci_wald[2] - ci_wald[1],
        ci_exact[2] - ci_exact[1],
        ci_tulap[2] - ci_tulap[1],
        ci_fima[2] - ci_fima[1]
      ),
      time_ms = c(time_wald, time_exact, time_tulap, time_fima)
    )
    
    pb$tick()  
  }
}

# Summarize results
summary_stats <- results %>%
  group_by(theta, method) %>%
  summarise(
    coverage = mean(coverage, na.rm = TRUE),
    avg_length = mean(ci_length, na.rm = TRUE),
    avg_time_ms = mean(time_ms, na.rm = TRUE),
    .groups = "drop"
  )

# View results
wald_stats <- summary_stats %>% filter(method == "Wald")
exact_stats <- summary_stats %>% filter(method == "Exact")
tulap_stats <- summary_stats %>% filter(method == "Tulap")
fima_stats <- summary_stats %>% filter(method == "FIMA")
