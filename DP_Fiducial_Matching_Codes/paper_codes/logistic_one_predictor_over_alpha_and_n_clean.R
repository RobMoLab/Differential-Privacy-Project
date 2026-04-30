################################################################################
# FIMA logistic regression simulation: one binary categorical predictor
# Grid: nominal level alpha x sample size n, with fixed privacy budget eps
################################################################################

# This script studies logistic regression with one binary predictor X1:
#   logit{P(Y = 1 | X1)} = beta0 + beta1 * X1.
# Since X1 is binary, beta0 and beta1 can be recovered from the two cell
# probabilities P(Y = 1 | X1 = 0) and P(Y = 1 | X1 = 1). FIMA is applied
# separately to each privatized cell proportion, and the resulting draws are
# transformed to draws of beta0 and beta1.

################################################################################
# Helper functions
################################################################################

# Numerically stable inverse-logit function.
inverse_logit <- function(eta) {
  1 / (1 + exp(-eta))
}

# Draw one Laplace(0, scale) random variable using inverse transform sampling.
# The input u is kept away from the exact endpoints to avoid log(0).
r_laplace <- function(scale) {
  u <- runif(1, min = -0.5, max = 0.5)
  u <- pmin(pmax(u, -0.5 + .Machine$double.eps), 0.5 - .Machine$double.eps)
  -scale * sign(u) * log(1 - 2 * abs(u))
}

# Privatize a sample proportion using additive Laplace noise with scale 1/(eps*n).
# If eps = Inf, no privacy noise is added.
privatize_proportion <- function(p_hat, eps, n) {
  if (is.infinite(eps)) return(p_hat)
  p_hat + r_laplace(scale = 1 / (eps * n))
}

# Generate H FIMA draws for one privatized proportion.
# If pi_priv = p_hat + Y, FIMA repeatedly samples an independent privacy noise Y*
# and forms delta = pi_priv - Y*, then samples from the fiducial Beta distribution.
fima_proportion <- function(pi_priv, H, eps, n, boundary_tol = 0.01, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  draws <- numeric(H)

  for (h in seq_len(H)) {
    delta_h <- if (is.infinite(eps)) {
      pi_priv
    } else {
      pi_priv - r_laplace(scale = 1 / (eps * n))
    }

    # Keep the beta-shape parameters strictly positive and avoid logits at 0 or 1.
    delta_h <- min(max(delta_h, boundary_tol), 1 - boundary_tol)
    draws[h] <- rbeta(1, n * delta_h + 0.5, n * (1 - delta_h) + 0.5)
  }

  draws
}

# Compute FIMA confidence intervals for beta0 and beta1 in the one-predictor model.
fima_logistic_one_predictor <- function(Y, X1, eps, H, alpha, seed = NULL) {
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  n0 <- sum(X1 == 0)
  n1 <- sum(X1 == 1)

  if (n0 == 0 || n1 == 0) {
    return(list(ci_beta0 = c(NA_real_, NA_real_), ci_beta1 = c(NA_real_, NA_real_)))
  }

  p0 <- mean(Y[X1 == 0])
  p1 <- mean(Y[X1 == 1])

  p0_priv <- privatize_proportion(p0, eps = eps, n = n0)
  p1_priv <- privatize_proportion(p1, eps = eps, n = n1)

  p0_fima <- fima_proportion(p0_priv, H = H, eps = eps, n = n0, seed = seed + 1)
  p1_fima <- fima_proportion(p1_priv, H = H, eps = eps, n = n1, seed = seed + 2)

  beta0_draws <- qlogis(p0_fima)
  beta1_draws <- qlogis(p1_fima) - beta0_draws

  list(
    ci_beta0 = quantile(beta0_draws, probs = c(alpha / 2, 1 - alpha / 2)),
    ci_beta1 = quantile(beta1_draws, probs = c(alpha / 2, 1 - alpha / 2))
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
    beta = c(beta0 = 0.5, beta1 = 2),
    p_x1 = 0.5,
    base_seed = 112
) {
  beta0_coverage <- matrix(NA_real_, length(sample_sizes), length(alpha_values),
                           dimnames = list(sample_sizes, alpha_values))
  beta1_coverage <- beta0_coverage
  beta0_CI_length <- beta0_coverage
  beta1_CI_length <- beta0_coverage

  for (a in seq_along(alpha_values)) {
    alpha <- alpha_values[a]

    for (m in seq_along(sample_sizes)) {
      n <- sample_sizes[m]

      # Fix the predictor design for all Monte Carlo replications at this n.
      set.seed(base_seed + m)
      X1 <- rbinom(n, size = 1, prob = p_x1)
      p_y <- inverse_logit(beta["beta0"] + beta["beta1"] * X1)

      coverage <- matrix(FALSE, B, length(beta))
      ci_length <- matrix(NA_real_, B, length(beta))

      for (b in seq_len(B)) {
        set.seed(b)
        Y <- rbinom(n, size = 1, prob = p_y)

        ci <- fima_logistic_one_predictor(
          Y = Y, X1 = X1, eps = eps, H = H, alpha = alpha,
          seed = base_seed + n + 10 * b
        )

        coverage[b, ] <- c(
          ci$ci_beta0[1] < beta["beta0"] && beta["beta0"] < ci$ci_beta0[2],
          ci$ci_beta1[1] < beta["beta1"] && beta["beta1"] < ci$ci_beta1[2]
        )
        ci_length[b, ] <- c(diff(ci$ci_beta0), diff(ci$ci_beta1))
      }

      beta0_coverage[m, a] <- mean(coverage[, 1], na.rm = TRUE)
      beta1_coverage[m, a] <- mean(coverage[, 2], na.rm = TRUE)
      beta0_CI_length[m, a] <- mean(ci_length[, 1], na.rm = TRUE)
      beta1_CI_length[m, a] <- mean(ci_length[, 2], na.rm = TRUE)
    }
  }

  list(
    beta0_coverage = beta0_coverage,
    beta1_coverage = beta1_coverage,
    beta0_CI_length = beta0_CI_length,
    beta1_CI_length = beta1_CI_length
  )
}

################################################################################
# Execute simulation
################################################################################

results <- run_simulation()
print(results)
