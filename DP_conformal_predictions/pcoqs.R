noisy_rc <- function(bounds, data, sigma) {
  
  # Noisy Range Count with Gaussian noise
  noisy_count <- sum(data >= bounds[1] & data <= bounds[2]) + rnorm(1, mean = 0, sd = sigma)
  
  return(max(0, floor(noisy_count)))  # Ensure non-negative integer count
  
}

priv_quant <- function(data, alpha, rho, lower_bound = 0, upper_bound = 1, delta = 1e-10) {
  
  # Differentially Private Quantile Approximation
  n <- length(data)
  L <- upper_bound - lower_bound
  sigma <- sqrt(ceiling(log2(L/delta)) / (2 * rho))
  m <- ceiling((1 - alpha) * (n + 1))
  
  left <- lower_bound
  right <- upper_bound
  
  N <- ceiling(log2(L/delta))
  
  for (i in 1:N) {
    
    mid <- (left + right) / 2
    c <- noisy_rc(c(lower_bound, mid), data, sigma)
    
    if (c < m) {
      
      left <- mid + delta
      
    } else {
      
      right <- mid
      
    }
    
  }
  
  return(round((left + right) / 2, 2))
  
}


set.seed(42)
D <- sort(runif(1000))  # Sorted dataset between 0 and 1
alpha <- 0.1
rho <- 1

priv_quant(D, alpha, rho)
quantile(D, probs = 1 - alpha)

pcoqs <- function(model, X_cal, Y_cal, X_test,
                          alpha = 0.1, rho = 1,
                          predict_fun = function(model, X, ...) {
                            predict(model, newdata = data.frame(X), ...)
                          },
                          score_fun = function(y, y_hat) abs(y - y_hat),
                          output_fun = function(y_hat, q_hat) {
                            list(prediction = y_hat, lower = y_hat - q_hat, upper = y_hat + q_hat)
                          },
                          lower_bound = 0, upper_bound = 1, delta = 1e-10,
                          predict_args = list()) {
  
  # Predict on calibration set
  Y_cal_hat <- do.call(predict_fun, c(list(model = model, X = X_cal), predict_args))
  
  # Compute nonconformity scores
  scores <- score_fun(Y_cal, Y_cal_hat)
  scores_sorted <- sort(scores)
  
  # Estimate DP quantile
  q_hat <- priv_quant(scores_sorted, alpha, rho, lower_bound, upper_bound, delta)
  
  # Predict on test set
  Y_test_hat <- do.call(predict_fun, c(list(model = model, X = X_test), predict_args))
  
  # Format output (interval, label, etc.)
  output_list <- output_fun(Y_test_hat, q_hat)
  output_df <- as.data.frame(output_list)
  
  return(list(
    output = output_df,
    quantile = q_hat
  ))
  
}




##############################################################################
# Regression
##############################################################################

# Load required package
set.seed(123)

# Simulate data
n <- 150
p <- 5
X <- matrix(rnorm(n * p), ncol = p)
colnames(X) <- paste0("X", 1:p)
beta <- runif(p, -2, 2)
Y <- X %*% beta + rnorm(n)

# Split into calibration and test
index_cal <- sample(1:n, ceiling(0.7*n))
index_test <- setdiff(1:n, index_cal)
X_cal <- X[index_cal, ]; Y_cal <- Y[index_cal]
X_test <- X[index_test, ]; Y_test <- Y[index_test]

# Fit linear model
df_train <- data.frame(Y = Y_cal, X_cal)
lm_model <- lm(Y ~ ., data = df_train)

# Run conformal prediction
result_reg <- pcoqs(
  model = lm_model,
  X_cal = X_cal,
  Y_cal = Y_cal,
  X_test = X_test,
  alpha = 0.1,
  rho = 0.1,
  score_fun = function(y, y_hat) abs(y - y_hat),  # absolute residual
  output_fun = function(y_hat, q_hat) {
    data.frame(
      prediction = y_hat,
      lower = y_hat - q_hat,
      upper = y_hat + q_hat
    )
  },
  lower_bound = 0,
  upper_bound = 10
)

# Evaluate coverage
intervals <- result_reg$output
covered <- (Y_test >= intervals$lower) & (Y_test <= intervals$upper)
cat("Regression coverage:", mean(covered), "\n")



##############################################################################
# Classification
##############################################################################

library(nnet)   # for multinom
library(dplyr)
set.seed(456)

# Simulate 3-class classification data
n <- 300; p <- 4
X <- matrix(rnorm(n * p), ncol = p)
colnames(X) <- paste0("X", 1:p)
class_probs <- t(apply(X[, 1:3], 1, function(row) {
  logits <- c(0, row[1] + row[2], row[3] - row[2])
  probs <- exp(logits) / sum(exp(logits))
  probs
}))
Y <- apply(class_probs, 1, function(p) sample(1:3, 1, prob = p))

# Split into calibration and test
index_cal <- sample(1:n, ceiling(0.7*n))
index_test <- setdiff(1:n, index_cal)
X_cal <- X[index_cal, ]; Y_cal <- Y[index_cal]
X_test <- X[index_test, ]; Y_test <- Y[index_test]

# Fit multinomial logistic regression
df_cal <- data.frame(Y = as.factor(Y_cal), X_cal)
model <- multinom(Y ~ ., data = df_cal, trace = FALSE)

# Predict softmax probabilities
predict_probs <- function(model, X, ...) {
  predict(model, newdata = data.frame(X), type = "probs")
}

# Nonconformity score: 1 - true class probability
prob_score <- function(y, probs) {
  sapply(1:length(y), function(i) 1 - probs[i, y[i]])
}

# Set-valued output: include all classes with prob >= 1 - q
set_output <- function(probs, q_hat) {
  sets <- apply(probs, 1, function(p) {
    labels <- which(p >= 1 - q_hat)
    paste(sort(labels), collapse = ",")
  })
  data.frame(predicted_set = sets)
}

# Run conformal prediction
result_mc <- pcoqs(
  model = model,
  X_cal = X_cal,
  Y_cal = Y_cal,
  X_test = X_test,
  alpha = 0.1,
  rho = 0.1,
  predict_fun = predict_probs,
  score_fun = prob_score,
  output_fun = set_output,
  lower_bound = 0,
  upper_bound = 1
)

# Evaluate set coverage
pred_sets <- strsplit(as.character(result_mc$output$predicted_set), ",")
covered <- mapply(function(set, true_y) {
  as.character(true_y) %in% set
}, pred_sets, Y_test)

cat("Multiclass coverage (set contains true label):", mean(covered), "\n")

