################################################################################
# FIMA logistic regression simulation: two binary categorical predictors
# Grid: nominal level alpha x sample size n, with fixed privacy budget eps
################################################################################

# This script studies logistic regression with two binary predictors:
#   logit{P(Y = 1 | X1, X2)} = beta0 + beta1 * X1 + beta2 * X2.
# For the additive model, beta0, beta1, and beta2 are obtained from the cell
# probabilities for (X1, X2) = (0,0), (1,0), and (0,1):
#   beta0 = logit(p00), beta1 = logit(p10) - beta0,
#   beta2 = logit(p01) - beta0.
# FIMA is applied independently to these three privatized cell proportions.

################################################################################
# Helper functions
################################################################################

inverse_logit <- function(eta) {
  1 / (1 + exp(-eta))
}

r_laplace <- function(scale) {
  u <- runif(1, min = -0.5, max = 0.5)
  u <- pmin(pmax(u, -0.5 + .Machine$double.eps), 0.5 - .Machine$double.eps)
  -scale * sign(u) * log(1 - 2 * abs(u))
}

privatize_proportion <- function(p_hat, eps, n) {
  if (is.infinite(eps)) return(p_hat)
  p_hat + r_laplace(scale = 1 / (eps * n))
}

fima_proportion <- function(pi_priv, H, eps, n, boundary_tol = 0.01, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  draws <- numeric(H)
  for (h in seq_len(H)) {
    delta_h <- if (is.infinite(eps)) {
      pi_priv
    } else {
      pi_priv - r_laplace(scale = 1 / (eps * n))
    }
    delta_h <- min(max(delta_h, boundary_tol), 1 - boundary_tol)
    draws[h] <- rbeta(1, n * delta_h + 0.5, n * (1 - delta_h) + 0.5)
  }
  draws
}

# Compute FIMA confidence intervals for beta0, beta1, and beta2.
fima_logistic_two_predictors <- function(Y, X1, X2, eps, H, alpha, seed = NULL) {
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  cell_00 <- X1 == 0 & X2 == 0
  cell_10 <- X1 == 1 & X2 == 0
  cell_01 <- X1 == 0 & X2 == 1

  n00 <- sum(cell_00)
  n10 <- sum(cell_10)
  n01 <- sum(cell_01)

  if (n00 == 0 || n10 == 0 || n01 == 0) {
    return(list(
      ci_beta0 = c(NA_real_, NA_real_),
      ci_beta1 = c(NA_real_, NA_real_),
      ci_beta2 = c(NA_real_, NA_real_)
    ))
  }

  p00 <- mean(Y[cell_00])
  p10 <- mean(Y[cell_10])
  p01 <- mean(Y[cell_01])

  p00_priv <- privatize_proportion(p00, eps = eps, n = n00)
  p10_priv <- privatize_proportion(p10, eps = eps, n = n10)
  p01_priv <- privatize_proportion(p01, eps = eps, n = n01)

  p00_fima <- fima_proportion(p00_priv, H = H, eps = eps, n = n00, seed = seed + 1)
  p10_fima <- fima_proportion(p10_priv, H = H, eps = eps, n = n10, seed = seed + 2)
  p01_fima <- fima_proportion(p01_priv, H = H, eps = eps, n = n01, seed = seed + 3)

  beta0_draws <- qlogis(p00_fima)
  beta1_draws <- qlogis(p10_fima) - beta0_draws
  beta2_draws <- qlogis(p01_fima) - beta0_draws

  list(
    ci_beta0 = quantile(beta0_draws, probs = c(alpha / 2, 1 - alpha / 2)),
    ci_beta1 = quantile(beta1_draws, probs = c(alpha / 2, 1 - alpha / 2)),
    ci_beta2 = quantile(beta2_draws, probs = c(alpha / 2, 1 - alpha / 2))
  )
}

################################################################################
# Main simulation
################################################################################

run_simulation <- function(
    alpha_values = c(0.01, 0.05, 0.10),
    sample_sizes = c(15, 30, 100, 200, 500, 1000, 2000),
    eps = 1,
    H = 10^4,
    B = 10^4,
    beta = c(beta0 = 0.5, beta1 = 2, beta2 = -2),
    p_x1 = 0.5,
    p_x2 = 0.5,
    base_seed = 144
) {
  beta0_coverage <- matrix(NA_real_, length(sample_sizes), length(alpha_values),
                           dimnames = list(sample_sizes, alpha_values))
  beta1_coverage <- beta0_coverage
  beta2_coverage <- beta0_coverage
  beta0_CI_length <- beta0_coverage
  beta1_CI_length <- beta0_coverage
  beta2_CI_length <- beta0_coverage

  for (a in seq_along(alpha_values)) {
    alpha <- alpha_values[a]

    for (m in seq_along(sample_sizes)) {
      n <- sample_sizes[m]

      # Fix the predictor design for all Monte Carlo replications at this n.
      set.seed(base_seed + m)
      X1 <- rbinom(n, size = 1, prob = p_x1)
      X2 <- rbinom(n, size = 1, prob = p_x2)
      p_y <- inverse_logit(beta["beta0"] + beta["beta1"] * X1 + beta["beta2"] * X2)

      coverage <- matrix(FALSE, B, length(beta))
      ci_length <- matrix(NA_real_, B, length(beta))

      for (b in seq_len(B)) {
        set.seed(b)
        Y <- rbinom(n, size = 1, prob = p_y)

        ci <- fima_logistic_two_predictors(
          Y = Y, X1 = X1, X2 = X2, eps = eps, H = H, alpha = alpha,
          seed = base_seed + n + 10 * b
        )

        coverage[b, ] <- c(
          ci$ci_beta0[1] < beta["beta0"] && beta["beta0"] < ci$ci_beta0[2],
          ci$ci_beta1[1] < beta["beta1"] && beta["beta1"] < ci$ci_beta1[2],
          ci$ci_beta2[1] < beta["beta2"] && beta["beta2"] < ci$ci_beta2[2]
        )
        ci_length[b, ] <- c(diff(ci$ci_beta0), diff(ci$ci_beta1), diff(ci$ci_beta2))
      }

      beta0_coverage[m, a] <- mean(coverage[, 1], na.rm = TRUE)
      beta1_coverage[m, a] <- mean(coverage[, 2], na.rm = TRUE)
      beta2_coverage[m, a] <- mean(coverage[, 3], na.rm = TRUE)
      beta0_CI_length[m, a] <- mean(ci_length[, 1], na.rm = TRUE)
      beta1_CI_length[m, a] <- mean(ci_length[, 2], na.rm = TRUE)
      beta2_CI_length[m, a] <- mean(ci_length[, 3], na.rm = TRUE)
    }
  }

  list(
    beta0_coverage = beta0_coverage,
    beta1_coverage = beta1_coverage,
    beta2_coverage = beta2_coverage,
    beta0_CI_length = beta0_CI_length,
    beta1_CI_length = beta1_CI_length,
    beta2_CI_length = beta2_CI_length
  )
}

################################################################################
# Execute simulation
################################################################################

results <- run_simulation()
print(results)
