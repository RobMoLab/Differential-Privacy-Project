#' FIMA Algorithm
#' 
#' @param pi0 Numeric: observed private proportion.
#' @param eps Numeric: Privacy budget parameter (must be > 0).
#' @param n Integer: Sample size (must be > 0).
#' @param H Integer: Number of samples to generate (default = 1).
#' @param seed Integer: Random seed for reproducibility (default = 123).
#'
#' @return Numeric vector of length `H` containing distribution of the estimator 
#'         

fima <- function(pi0, eps, n, H = 1, delta = 1, seed = 123) {
  
  set.seed(seed)  # Ensure reproducibility
  
  # Step 1: Generate random weights uniformly in [-0.5, 0.5]
  Wj <- runif(H, -0.5, 0.5)
  
  # Step 2: "denoise" pivatized proportion
  pi_star <- pi0 + (delta / (eps * n)) * sign(Wj) * log(1 - 2 * abs(Wj))
  
  # Initialize output vector
  theta_ifb <- rep(NA, H)
  
  # Step 3: Sample from Beta distribution 
  valid_samples <- (pi_star < 1) & (pi_star > 0)
  theta_ifb[valid_samples] <- rbeta(
    sum(valid_samples),
    shape1 = n * pi_star[valid_samples] + 0.5,
    shape2 = n * (1 - pi_star[valid_samples]) + 0.5
  )
  
  # Handle edge cases (pi_star <= 0 or >= 1)
  theta_ifb[pi_star >= 1] <- 1 - .Machine$double.eps  # Clip to just below 1
  theta_ifb[pi_star <= 0] <- .Machine$double.eps      # Clip to just above 0
  
  return(theta_ifb)
}

# Example
fima(pi0 =0.3, eps = 1, H=10, n=100)
