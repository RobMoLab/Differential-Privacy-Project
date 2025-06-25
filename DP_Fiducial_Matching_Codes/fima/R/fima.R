#' Produce DP proportion with additive Laplace mechanism
#' @description This function is a wrapper to directly extract p-values (where possible) and confidence intervals (where needed) from the FIMA distributions produced by functions \code{fima_prop}, \code{fima_2prop}, \code{fima_chi2} and \code{fima_logit} (see Romanus et al., 2025).
#' @param obj The object output from the functions \code{fima_prop}, \code{fima_2prop}, \code{fima_chi2} and \code{fima_logit} (it consists in a \code{vector} representing the FIMA distribution).
#' @param theta0 A \code{double} value representing the parameter value to test under the null hypothesis (default value is \code{theta0 = 0.5}). This is used only when \code{obj} comes from \code{fima_prop} or \code{fima_2prop}.
#' @param ha A \code{character} representing the one-sided alternative for the hypothesis test. The options are \code{"greater"} (i.e. we want to test if the true parameter is greater than \code{theta0}) or \code{"smaller"} (i.e. we want to test if the true parameter is smaller than \code{theta0}).
#' @param alpha A \code{double} value representing the significance level for the \code{1 - alpha} confidence intervals on the FIMA distribution (default value is \code{alpha = 0.05}).
#' @return A \code{list} containing the following objects:
#' \describe{
#'  \item{null}{The value of the parameter tested under the null hypothesis. This is produced for objects output from \code{fima_prop} and \code{fima_2prop} functions.}
#'  \item{alternative}{The direction of the hypothesis test under the alternative. This is produced for objects output from \code{fima_prop} and \code{fima_2prop} functions.}
#'  \item{p_value}{The p-value of the hypothesis test. This is produced for objects output from \code{fima_prop}, \code{fima_2prop} and \code{fima_chi2} functions.}
#'  \item{conf_int}{A \code{double} vector containing the lower and upper bound of the confidence intervals for the parameter of interest. This is produced for objects output from \code{fima_prop}, \code{fima_2prop} and \code{fima_logit} functions.}
#' }
#' @author Roberto Molinari and Ogonnaya M. Romanus
#' @import stats
#' @import utils
#' @export
#' @examples
#' \dontrun{
#' # Inference for one-sample proportion
#' p <- 0.8 # true proportion
#' n <- 30 # sample size
#' set.seed(14) # seed for reproducibility
#' x <- rbinom(n, 1, prob = p) # simulate data
#' eps <- 1 # epsilon-DP privacy budget
#' pi <- dp_prop(mean(x), eps = eps, n = n) # produce DP proportion
#' H <- 10^4 # number of simulations for FIMA distribution
#' dist <- fima_prop(pi, eps = eps, n = n, H = H) # produce FIMA distribution
#' fima_infer(dist) # obtain inferential quantities
#' }
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

fima_prop <- function(dp_stat, n, eps = 1, delta = 1, H = 10^4, seed = NA) {

  if(!is.na(seed)) set.seed(seed)

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

fima_count <- function(dp_stat, n, eps = 1, delta = 1, H = 10^4, terms = 1, seed = NA) {

  if(!is.na(seed)) set.seed(seed)

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

fima_2prop <- function(dp_stat1, dp_stat2, n1, n2, eps = 1, delta = 1, H = 10^4, seed = NA) {

  # Generate distributions for both groups
  fima_prop1 <- fima_prop(dp_stat = dp_stat1, n = n1, eps = eps, delta = delta, H = H, seed = seed)
  fima_prop2 <- fima_prop(dp_stat = dp_stat2, n = n2, eps = eps, delta = delta, H = H, seed = seed)

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

fima_chi2 <- function(dp_table, n, eps = 1, delta = 2, H = 10^4, seed = NA) {

  # Compute chi-2 statistic on DP table
  dp_table <- pmax(dp_table, 0) + 0.5 # Haldane-Anscombe Correction
  n_star <- sum(dp_table)
  marginal_row <- apply(dp_table, 1, sum)
  marginal_col <- apply(dp_table, 2, sum)
  expected <- outer(marginal_row/n_star, marginal_col/n_star)*n_star
  chi2_obs <- sum((dp_table - expected)^2 / expected)

  # Vectorized FIMA for marginals
  fima_marginal_row <- sapply(marginal_row, fima_count, n = n, eps = eps, delta = delta, H = H, terms = nrow(dp_table), seed = seed)
  fima_marginal_col <- sapply(marginal_col, fima_count, n = n, eps = eps, delta = delta, H = H, terms = ncol(dp_table), seed = seed)

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


#' FIMA Inference
#' @description This function is a wrapper to directly extract p-values (where possible) and confidence intervals (where needed) from the FIMA distributions produced by functions \code{fima_prop}, \code{fima_2prop}, \code{fima_chi2} and \code{fima_logit} (see Romanus et al., 2025).
#' @param obj The object output from the functions \code{fima_prop}, \code{fima_2prop}, \code{fima_chi2} and \code{fima_logit} (it consists in a \code{vector} representing the FIMA distribution).
#' @param theta0 A \code{double} value representing the parameter value to test under the null hypothesis (default value is \code{theta0 = 0.5}). This is used only when \code{obj} comes from \code{fima_prop} or \code{fima_2prop}.
#' @param ha A \code{character} representing the one-sided alternative for the hypothesis test. The options are \code{"greater"} (i.e. we want to test if the true parameter is greater than \code{theta0}) or \code{"smaller"} (i.e. we want to test if the true parameter is smaller than \code{theta0}).
#' @param alpha A \code{double} value representing the significance level for the \code{1 - alpha} confidence intervals on the FIMA distribution (default value is \code{alpha = 0.05}).
#' @return A \code{list} containing the following objects:
#' \describe{
#'  \item{null}{The value of the parameter tested under the null hypothesis. This is produced for objects output from \code{fima_prop} and \code{fima_2prop} functions.}
#'  \item{alternative}{The direction of the hypothesis test under the alternative. This is produced for objects output from \code{fima_prop} and \code{fima_2prop} functions.}
#'  \item{p_value}{The p-value of the hypothesis test. This is produced for objects output from \code{fima_prop}, \code{fima_2prop} and \code{fima_chi2} functions.}
#'  \item{conf_int}{A \code{double} vector containing the lower and upper bound of the confidence intervals for the parameter of interest. This is produced for objects output from \code{fima_prop}, \code{fima_2prop} and \code{fima_logit} functions.}
#' }
#' @author Roberto Molinari and Ogonnaya M. Romanus
#' @import stats
#' @import utils
#' @export
#' @examples
#' \dontrun{
#' # Inference for one-sample proportion
#' p <- 0.8 # true proportion
#' n <- 30 # sample size
#' set.seed(14) # seed for reproducibility
#' x <- rbinom(n, 1, prob = p) # simulate data
#' eps <- 1 # epsilon-DP privacy budget
#' pi <- dp_prop(mean(x), eps = eps, n = n) # produce DP proportion
#' H <- 10^4 # number of simulations for FIMA distribution
#' dist <- fima_prop(pi, eps = eps, n = n, H = H) # produce FIMA distribution
#' fima_infer(dist) # obtain inferential quantities
#' }
fima_infer <- function(obj, theta0 = 0.5, ha = "greater", alpha = 0.05) {

  if(class(obj) == "prop") {

    if(ha == "greater") {

      p_val <- (sum(obj > theta0) + 1)/(length(obj) + 1)

    } else {

      p_val <- (sum(obj < theta0) + 1)/(length(obj) + 1)

    }

    ci <- quantile(obj, probs = c(alpha/2, 1 - alpha/2))

    return(list("null" = theta0, "alternative" = ha, "p_value" = p_val, "conf_int" = ci))

  } else if(class(obj) == "2prop") {

    if(ha == "smaller") {

      p_val <- (sum(obj > theta0) + 1)/(length(obj) + 1)

    } else {

      p_val <- (sum(obj < theta0) + 1)/(length(obj) + 1)

    }

    ci <- quantile(obj, probs = c(alpha/2, 1 - alpha/2))

    return(list("null" = theta0, "alternative" = ha, "p_value" = p_val, "conf_int" = ci))

  } else if (class(obj) == "chi2") {

    p_val <- (sum(obj$dist >= obj$obs) + 1)/(length(obj$dist) + 1)

    return(list("p_value" = p_val))

  } else if(class(obj) == "logistic"){

    ci <- apply(obj, 2, quantile, probs = c(alpha/2, 1 - alpha/2))

    return(list("conf_int" = ci))

  }

}
