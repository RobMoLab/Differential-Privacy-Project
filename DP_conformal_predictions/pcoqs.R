NoisyRC <- function(range_bounds, D, sigma) {
  # Noisy Range Count with Gaussian noise
  a <- range_bounds[1]
  b <- range_bounds[2]
  
  count <- sum(D >= a & D <= b)
  noise <- rnorm(1, mean = 0, sd = sigma)
  noisy_count <- count + noise
  
  return(max(0, floor(noisy_count)))  # Ensure non-negative integer count
}

PrivQuant <- function(D, alpha, rho, lower_bound = 0, upper_bound = 1, delta = 1e-10) {
  # Differentially Private Quantile Approximation
  
  n <- length(D)
  L <- upper_bound - lower_bound
  sigma <- sqrt(ceiling(log2(L)) / (2 * rho))
  m <- ceiling((1 - alpha) * (n + 1))
  
  left <- lower_bound
  right <- upper_bound
  
  N <- ceiling(log2(L/delta))
  
  for (i in 1:N) {
    
    mid <- (left + right) / 2
    c <- NoisyRC(c(lower_bound, mid), D, sigma)
    
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
alpha <- 0.25
rho <- 1.0

PrivQuant(D, alpha, rho)
quantile(D, probs = 1 - alpha)


pcoq <- function(X_train, Y_train,
                                 X_cal, Y_cal,
                                 X_test,
                                 alpha = 0.1,
                                 rho = 1.0,
                                 model_fun = function(X, Y) { lm(Y ~ ., data = data.frame(Y, X)) },
                                 predict_fun = function(model, X) { predict(model, newdata = data.frame(X)) },
                                 lower_bound = 0, upper_bound = 1, delta = 1e-10) {
  
  # Fit model on training data
  model <- model_fun(X_train, Y_train)
  
  # Predict on calibration set
  Y_cal_hat <- predict_fun(model, X_cal)
  
  # Compute conformal scores
  scores <- abs(Y_cal - Y_cal_hat)
  scores_sorted <- sort(scores)  # required by PrivQuant
  
  # Estimate DP quantile
  q_hat <- PrivQuant(scores_sorted, alpha, rho, lower_bound, upper_bound, delta)
  
  # Predict on test set
  Y_test_hat <- predict_fun(model, X_test)
  
  # Return predictive intervals
  intervals <- data.frame(
    lower = Y_test_hat - q_hat,
    upper = Y_test_hat + q_hat,
    prediction = Y_test_hat
  )
  
  return(list(
    intervals = intervals,
    quantile = q_hat,
    model = model
  ))
}



# 1. Generate synthetic data
set.seed(42)
n <- 100
p <- 5

X <- matrix(rnorm(n * p), nrow = n)
beta <- runif(p, -2, 2)
Y <- X %*% beta + rnorm(n)

# 2. Split into train / calibration / test
train_idx <- 1:60
cal_idx   <- 61:90
test_idx  <- 91:100

X_train <- X[train_idx, ]
Y_train <- Y[train_idx]

X_cal <- X[cal_idx, ]
Y_cal <- Y[cal_idx]

X_test <- X[test_idx, ]

# 3. Run DP conformal prediction
result <- pcoq(X_train, Y_train, X_cal, Y_cal, X_test,
                               alpha = 0.1, rho = 1.0,
                               lower_bound = 0, upper_bound = 10, delta = 1e-5)

# 4. View predictive intervals
head(result$intervals)


# Extract test set true responses
Y_test <- Y[test_idx]

# Extract prediction intervals
intervals <- result$intervals

# Check coverage: is Y_test inside [lower, upper]?
coverage <- (Y_test >= intervals$lower) & (Y_test <= intervals$upper)

# Print overall coverage rate
cat("Coverage rate:", mean(coverage), "\n")

# Optional: View which points are covered
data.frame(
  Y_test = round(Y_test, 2),
  Lower = round(intervals$lower, 2),
  Upper = round(intervals$upper, 2),
  Covered = coverage
)
