#' Differentially Private Proportion with Additive Laplace Mechanism
#' @description Based on a sample proportion, the function outputs a differentially private proportion by adding calibrated Laplace noise. More specifically, this guarantees pure DP (with privacy budget \eqn{\epsilon})
#' @param prop A \code{double} value (scalar or vector) representing the sample proportion(s) that need(s) to be privatized.
#' @param eps A \code{double} value representing the privacy budget \eqn{\epsilon} to guarantee \eqn{\epsilon}-DP for the resulting proportion (default value is \code{eps = 1}).
#' @param delta A \code{double} representing the sensitivity of the count associated to the proportion. As a consequence the variance of the Laplace noise is \eqn{\Delta/(n \epsilon)} (default value is \code{delta = 1}).
#' @return A \code{double} (scalar or vector) representing the privatized proportion(s) guaranteeing \eqn{\epsilon}-DP (for each element in case \code{prop} is a vector).
#' @author Roberto Molinari and Ogonnaya M. Romanus
#' @import stats
#' @import utils
#' @export
#' @examples
#' \dontrun{
#' # Differentially private proportion
#' p <- 0.8 # true probability
#' n <- 30 # sample size
#' set.seed(14) # seed for reproducibility
#' x <- rbinom(n, 1, prob = p) # simulate data
#' eps <- 1 # epsilon-DP privacy budget
#' pi <- dp_prop(mean(x), eps = eps, n = n) # produce DP proportion
#' }
dp_prop <- function(prop, eps = 1, n, delta = 1) {

  W <- runif(length(prop), -0.5, 0.5)
  dp_stat <- prop + (delta / (eps * n)) * sign(W) * log(1 - 2 * abs(W))

  return(dp_stat)

}

#' Differentially Private Count with Additive Laplace Mechanism
#' @description Based on a sample count, the function outputs a differentially private count by adding calibrated Laplace noise. More specifically, this guarantees pure DP (with privacy budget \eqn{\epsilon})
#' @param count A \code{double} value (scalar or vector) representing the sample count(s) that need(s) to be privatized.
#' @param eps A \code{double} value representing the privacy budget \eqn{\epsilon} to guarantee \eqn{\epsilon}-DP for the resulting proportion (default value is \code{eps = 1}).
#' @param delta A \code{double} representing the sensitivity of the count. As a consequence the variance of the Laplace noise is \eqn{\Delta/(\epsilon)} (default value is \code{delta = 1}).
#' @return A \code{double} (scalar or vector) representing the privatized count(s) guaranteeing \eqn{\epsilon}-DP (for each element in case \code{count} is a vector).
#' @author Roberto Molinari and Ogonnaya M. Romanus
#' @import stats
#' @import utils
#' @export
#' @examples
#' \dontrun{
#' # Differentially private 2x2 table of independent counts
#' joint_probs <- c(0.25, 0.25, 0.25, 0.25) # true probabilities for multinomial
#' n <- 100 # total sample size
#' set.seed(124) # seed for reproducibility
#' counts <- rmultinom(1, n, prob = joint_probs) # sample counts for each category
#' eps <- 1 # epsilon-DP privacy budget
#' priv_counts <- dp_count(counts, eps, delta = 2) # DP counts for each category
#' priv_tab <- matrix(priv_counts, 2, 2) # Contingency table of DP counts
#' }
dp_count <- function(count, eps = 1, delta = 1) {

  W <- runif(length(count), -0.5, 0.5)
  dp_noise <- (delta / eps) * sign(W) * log(1 - 2 * abs(W))
  dp_stat <- count + dp_noise

  return(dp_stat)

}


#' FIMA for One-Sample Proportion
#' @description Provides the FIMA distribution for inference on the true success probability in the one-sample proportion setting (see Romanus et al., 2025).
#' @param dp_stat The \code{double} value representing a privatized proportion with Laplace noise.
#' @param n A \code{double} value representing the sample size on which the non-private proportion was computed.
#' @param eps A \code{double} value representing the privacy budget \eqn{\epsilon} used to guarantee \eqn{\epsilon}-DP for the proportion in \code{dp_stat} (default value is \code{eps = 1}).
#' @param delta A \code{double} representing the count-sensitivity used to guarantee \eqn{\epsilon}-DP for the proportion in \code{dp_stat} (default value is \code{delta = 1}).
#' @param H An \code{integer} representing the number of bootstrap solutions to approximate the FIMA distribution (default value is \code{H = 10^4}).
#' @param seed An \code{integer} representing a seed value to generate the random variables used in the FIMA to generate solutions (default value is \code{seed = NA}, hence no seed).
#' @return A \code{vector} containing the sequence of FIMA solutions to be used for inference on the true success probability
#' @author Roberto Molinari and Ogonnaya M. Romanus
#' @import stats
#' @import utils
#' @export
#' @examples
#' \dontrun{
#' # FIMA distribution for one-sample proportion
#' p <- 0.8 # true proportion
#' n <- 30 # sample size
#' set.seed(14) # seed for reproducibility
#' x <- rbinom(n, 1, prob = p) # simulate data
#' eps <- 1 # epsilon-DP privacy budget
#' pi <- dp_prop(mean(x), eps = eps, n = n) # produce DP proportion
#' H <- 10^4 # number of simulations for FIMA distribution
#' dist <- fima_prop(pi, eps = eps, n = n, H = H) # produce FIMA distribution
#' }
#' @references Romanus, O.M., Boulaguiem, Y. & Molinari, R. (2025) "Fiducial Matching: Differentially Private Inference for Categorical Data".
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

#' FIMA for One-Sample Count
#' @description Provides the FIMA distribution for inference on the true success probability in the one-sample count setting (see Romanus et al., 2025).
#' @param dp_stat The \code{double} value representing a privatized count with Laplace noise. This is mainly used within the \code{fima_chi2} function.
#' @param n A \code{double} value representing the sample size on which the non-private proportion was computed.
#' @param eps A \code{double} value representing the privacy budget \eqn{\epsilon} used to guarantee \eqn{\epsilon}-DP for the proportion in \code{dp_stat} (default value is \code{eps = 1}).
#' @param delta A \code{double} representing the count-sensitivity used to guarantee \eqn{\epsilon}-DP for the proportion in \code{dp_stat} (default value is \code{delta = 1}).
#' @param H An \code{integer} representing the number of bootstrap solutions to approximate the FIMA distribution (default value is \code{H = 10^4}).
#' @param terms An \code{integer} representing the number of DP counts summed to obtain the final DP count \code{dp_stat} (default value is \code{terms = 1}, hence the DP count \code{dp_stat} is a direct output of a single private count).
#' @param seed An \code{integer} representing a seed value to generate the random variables used in the FIMA to generate solutions (default value is \code{seed = NA}, hence no seed).
#' @return A \code{vector} containing the sequence of FIMA solutions to be used for inference on the true success probability
#' @author Roberto Molinari and Ogonnaya M. Romanus
#' @import stats
#' @import utils
#' @export
#' @examples
#' \dontrun{
#' # FIMA distribution for one-sample count setting
#' p <- 0.8 # true proportion
#' n <- 30 # sample size
#' set.seed(14) # seed for reproducibility
#' x <- rbinom(n, 1, prob = p) # simulate data
#' eps <- 1 # epsilon-DP privacy budget
#' pi <- dp_count(x, eps = eps, n = n) # produce DP count
#' H <- 10^4 # number of simulations for FIMA distribution
#' dist <- fima_count(pi, eps = eps, n = n, H = H) # produce FIMA distribution for probability
#' }
#' @references Romanus, O.M., Boulaguiem, Y. & Molinari, R. (2025) "Fiducial Matching: Differentially Private Inference for Categorical Data".
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

#' FIMA for Independent Two-Sample Proportions
#' @description Provides the FIMA distribution for inference on the difference of true success probabilities in the independent two-sample proportions setting (see Romanus et al., 2025).
#' @param dp_stat1 The \code{double} value representing a privatized proportion with Laplace noise for the first sample.
#' @param dp_stat2 The \code{double} value representing a privatized proportion with Laplace noise for the second sample.
#' @param n1 An \code{integer} value representing the size of the first sample on which the proportion was computed.
#' @param n2 An \code{integer} value representing the size of the second sample on which the proportion was computed.
#' @param eps A \code{double} value representing the privacy budget \eqn{\epsilon} used to guarantee \eqn{\epsilon}-DP for the proportions in \code{dp_stat1} and \code{dp_stat2} (default value is \code{eps = 1}).
#' @param delta A \code{double} representing the count-sensitivity used to guarantee \eqn{\epsilon}-DP for the proportion in \code{dp_stat} (default value is \code{delta = 1}).
#' @param H An \code{integer} representing the number of bootstrap solutions to approximate the FIMA distribution (default value is \code{H = 10^4}).
#' @param seed An \code{integer} representing a seed value to generate the random variables used in the FIMA to generate solutions (default value is \code{seed = NA}, hence no seed).
#' @return A \code{vector} containing the sequence of FIMA solutions to be used for inference on the true difference in success probabilities
#' @author Roberto Molinari and Ogonnaya M. Romanus
#' @import stats
#' @import utils
#' @export
#' @examples
#' \dontrun{
#' # FIMA distribution for two-sample proportions
#' p1 <- 0.8 # true success probability for first sample
#' p2 <- 0.9 # true success probability for second sample
#' n <- 100 # sample size (common for this example)
#' set.seed(14) # seed for reproducibility
#' x1 <- rbinom(n, 1, prob = p1) # simulate data for first sample
#' x2 <- rbinom(n, 1, prob = p2) # simulate data for second sample
#' pi1 <- dp_prop(mean(x1), eps = 1, n = n) # produce DP proportion for first sample
#' pi2 <- dp_prop(mean(x2), eps = 1, n = n) # produce DP proportion for second sample
#' dist <- fima_2prop(pi1, pi2, n1 = n, n2 = n, eps = eps, H = H) # produce FIMA distribution for difference in probabilities
#' }
#' @references Romanus, O.M., Boulaguiem, Y. & Molinari, R. (2025) "Fiducial Matching: Differentially Private Inference for Categorical Data".
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


#' FIMA for \eqn{\chi^2}-test
#' @description Provides the FIMA distribution of the \eqn{\chi^2} statistic under the null hypothesis of independence between categorical variables (see Romanus et al., 2025).
#' @param dp_table A \code{double} matrix representing the contingency table containing the privatized counts with Laplace noise.
#' @param n An \code{integer} value representing the sample size (i.e. total original counts).
#' @param eps A \code{double} value representing the privacy budget \eqn{\epsilon} used to guarantee \eqn{\epsilon}-DP for the counts in \code{dp_table} (default value is \code{eps = 1}).
#' @param delta A \code{double} representing the count-sensitivity used to guarantee \eqn{\epsilon}-DP for the counts in \code{dp_table} (default value is \code{delta = 2}).
#' @param H An \code{integer} representing the number of bootstrap solutions to approximate the FIMA distribution (default value is \code{H = 10^4}).
#' @param seed An \code{integer} representing a seed value to generate the random variables used in the FIMA to generate solutions (default value is \code{seed = NA}, hence no seed).
#' @return A \code{vector} containing the sequence of FIMA solutions for the \eqn{\chi^2}-statistic under the null hypothesis of independence
#' @author Roberto Molinari and Ogonnaya M. Romanus
#' @import stats
#' @import utils
#' @export
#' @examples
#' \dontrun{
#' # FIMA distribution for two-sample proportions
#' joint_probs <- c(0.25, 0.25, 0.25, 0.25) # true probabilities for multinomial
#' n <- 100 # total sample size
#' set.seed(124) # seed for reproducibility
#' counts <- rmultinom(1, n, prob = joint_probs) # sample counts for each category
#' eps <- 1 # epsilon-DP privacy budget
#' priv_counts <- dp_count(counts, eps, delta = 2) # DP counts for each category
#' priv_tab <- matrix(priv_counts, 2, 2) # Contingency table of DP counts
#' dist <- fima_chi2(priv_tab, n, eps) # produce FIMA distribution for chi-2 statistic under independence
#' }
#' @references Romanus, O.M., Boulaguiem, Y. & Molinari, R. (2025) "Fiducial Matching: Differentially Private Inference for Categorical Data".
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
#' @description A wrapper function to extract p-values (where possible) and confidence intervals (where needed) from the FIMA distributions produced by functions \code{fima_prop}, \code{fima_2prop}, \code{fima_chi2} and \code{fima_logit} (see Romanus et al., 2025).
#' @param obj The object output from the functions \code{fima_prop}, \code{fima_2prop}, \code{fima_chi2} and \code{fima_logit} (it consists in a \code{vector} representing the FIMA distribution).
#' @param theta0 A \code{double} value representing the parameter value to test under the null hypothesis (default value is \code{theta0 = 0.5}). This is used only when \code{obj} comes from \code{fima_prop} or \code{fima_2prop}.
#' @param ha A \code{character} representing the one-sided alternative for the hypothesis test. The options are \code{"greater"} (i.e. we want to test if the true parameter is greater than \code{theta0}) or \code{"smaller"} (i.e. we want to test if the true parameter is smaller than \code{theta0}).
#' @param alpha A \code{double} value representing the significance level for the \eqn{1 - \alpha} confidence intervals on the FIMA distribution (default value is \code{alpha = 0.05}).
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
#' @references Romanus, O.M., Boulaguiem, Y. & Molinari, R. (2025) "Fiducial Matching: Differentially Private Inference for Categorical Data".
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
