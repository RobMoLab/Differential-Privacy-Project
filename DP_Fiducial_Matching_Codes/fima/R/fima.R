dp_prop <- function(true_prop, eps = 1, n, delta = 1) {

  W <- runif(length(true_prop), -0.5, 0.5)
  dp_stat <- true_prop + (delta / (eps * n)) * sign(W) * log(1 - 2 * abs(W))

  return(dp_stat)

}

dp_count <- function(count, eps = 1, delta = 1) {

  W <- runif(length(count), -0.5, 0.5)
  dp_noise <- (delta / eps) * sign(W) * log(1 - 2 * abs(W))
  dp_stat <- count + dp_noise

  return(dp_stat)

}

fima_prop <- function(dp_stat, n, eps = 1, delta = 1, H = 10^4, seed = 123) {

  set.seed(seed)

  # Seeds for Laplace privacy noise
  W <- runif(H, -0.5, 0.5)

  # Quantity on which to inverse binomial CDF
  pi_star <- dp_stat - (delta / (eps * n)) * sign(W) * log(1 - 2 * abs(W))

  # Initialize output vector
  theta <- rep(NA, H)

  # Sample from Beta distribution
  valid_samples <- (pi_star < 1) & (pi_star > 0) # Discard edge cases
  theta[valid_samples] <- rbeta(
    sum(valid_samples),
    shape1 = n * pi_star[valid_samples] + 0.5,
    shape2 = n * (1 - pi_star[valid_samples]) + 0.5
  )

  # Handle edge cases (pi_star <= 0 or >= 1)
  theta[pi_star >= 1] <- 1 - .Machine$double.eps  # Clip to just below 1
  theta[pi_star <= 0] <- .Machine$double.eps      # Clip to just above 0

  class(theta) <- "prop"

  return(theta)

}

fima_count <- function(dp_stat, n, eps = 1, delta = 1, H = 10^4, terms = 1, seed = 123) {

  set.seed(seed)

  # Seeds for Laplace privacy noise
  W <- runif(H*terms, -0.5, 0.5)
  vec_noise <- (delta / eps) * sign(W) * log(1 - 2 * abs(W))
  vec_noise <- sapply(split(vec_noise, ceiling(seq_along(vec_noise) / terms)), sum)

  # Quantity on which to inverse binomial CDF
  pi_star <- (dp_stat - vec_noise) / n

  # Initialize output vector
  theta <- rep(NA, H)

  # Sample from Beta distribution
  valid_samples <- (pi_star < 1) & (pi_star > 0) # Discard edge cases
  theta[valid_samples] <- rbeta(
    sum(valid_samples),
    shape1 = n * pi_star[valid_samples] + 0.5,
    shape2 = n * (1 - pi_star[valid_samples]) + 0.5
  )

  # Handle edge cases (pi_star <= 0 or >= 1)
  theta[pi_star >= 1] <- 1 - .Machine$double.eps  # Clip to just below 1
  theta[pi_star <= 0] <- .Machine$double.eps      # Clip to just above 0

  return(theta)

}

fima_2prop <- function(dp_stat1, dp_stat2, n1, n2, eps = 1, delta = 1, H = 10^4, seed = 123) {

  # Generate distributions for both groups
  fima_prop1 <- fima_prop(dp_stat = dp_stat1, n = n1, eps = eps, delta = delta, H = H, seed = seed + 1)
  fima_prop2 <- fima_prop(dp_stat = dp_stat2, n = n2, eps = eps, delta = delta, H = H, seed = seed + 2)

  # Compute differences
  theta <- fima_prop1 - fima_prop2

  class(theta) <- "2prop"

  return(theta)

}

chi2 <- function(tab){

  n <- sum(tab)
  marginal_row <- apply(tab, 1, sum)/n
  maringal_col <- apply(tab, 2, sum)/n
  expected <- outer(marginal_row, maringal_col) * n

  return(sum((tab - expected)^2 / expected))

}

fima_chi2 <- function(dp_table, n, eps = 1, delta = 2, H = 10^4, seed = 123) {

  # Compute chi-2 statistic on DP table
  dp_table <- pmax(dp_table, 0) + 0.5 # Haldane-Anscombe Correction
  n_star <- sum(dp_table)
  marginal_row <- apply(dp_table, 1, sum)
  marginal_col <- apply(dp_table, 2, sum)
  expected <- outer(marginal_row/n_star, marginal_col/n_star)*n_star
  chi2_obs <- sum((dp_table - expected)^2 / expected)

  # Vectorized FIMA for marginals
  fima_marginal_row <- sapply(marginal_row, fima_count, n = n, eps = eps, delta = delta, H = H, terms = nrow(dp_table), seed = seed + 1)
  fima_marginal_col <- sapply(marginal_col, fima_count, n = n, eps = eps, delta = delta, H = H, terms = ncol(dp_table), seed = seed + 2)

  # Generate FIMA joint probabilities and counts
  fima_counts <- array(NA, dim = c(H, nrow(dp_table), ncol(dp_table)))
  for (h in 1:H) {

    fima_joint_probs <- outer(fima_marginal_row[h, ], fima_marginal_col[h, ])
    fima_counts[h, , ] <- matrix(rmultinom(n = 1, size = n, prob = fima_joint_probs), nrow(dp_table), ncol(dp_table)) + 0.5 # Haldane-Anscombe Correction

  }

  fima_chi2_dist <- apply(dp_count(fima_counts, eps = eps, delta = delta), 1, chi2)

  out <- list("dist" = fima_chi2_dist, "obs" = chi2_obs)

  class(out) <- "chi2"

  return(out)

}

logit <- function(theta) {

  return(log(theta / (1 - theta)))

}

expit <- function(x) {

  return(exp(x) / (1 + exp(x)))

}

### This is wrong since each proportion has different base sample size n: needs to be fixed
fima_logistic <- function(dp_pi, n, eps = 1, delta = 2, H = 10^4, seed = 123) {

  # Generate distributions for betas
  fima_beta0 <- logit(fima_prop(dp_stat = dp_pi[1], n = n, eps = eps / length(dp_pi), delta = delta, H = H, seed = seed + 1))
  fima_betas <- logit(sapply(dp_pi[-1], fima_prop, n = n, eps = eps / length(dp_pi), delta = delta, H = H, seed = seed + 2)) - fima_beta0

  # Compute differences
  beta <- cbind(fima_beta0, fima_betas)
  colnames(beta) <- c("Beta_0", paste0("Beta_", 1:ncol(fima_betas)))

  class(beta) <- "logistic"

  return(beta)

}

fima_infer <- function(obj, theta0 = 0.5, ha = "greater", alpha = 0.05) {

  if(class(obj) == "prop") {

    if(ha == "greater") {

      p_val <- (sum(obj > theta0) + 1)/(length(obj) + 1)

    } else {

      p_val <- (sum(obj < theta0) + 1)/(length(obj) + 1)

    }

    ci <- quantile(obj, probs = c(alpha/2, 1 - alpha/2))

    return(list("null" = theta0, "alternative" = ha, "p-value" = p_val, "conf_int" = ci))

  } else if(class(obj) == "2prop") {

    if(ha == "smaller") {

      p_val <- (sum(obj > theta0) + 1)/(length(obj) + 1)

    } else {

      p_val <- (sum(obj < theta0) + 1)/(length(obj) + 1)

    }

    ci <- quantile(obj, probs = c(alpha/2, 1 - alpha/2))

    return(list("null" = theta0, "alternative" = ha, "p-value" = p_val, "conf_int" = ci))

  } else if (class(obj) == "chi2") {

    p_val <- (sum(obj$dist >= obj$obs) + 1)/(H + 1)

    return(list("p-value" = p_val))

  } else if(class(obj) == "logistic"){

    ci <- apply(obj, 2, quantile, probs = c(alpha/2, 1 - alpha/2))

    return(list("conf_int" = ci))

  }

}
