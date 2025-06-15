
# Chi-squared
diff <- 0.01
joint_probs0 <- c(0.25, 0.25, 0.25, 0.25) + diff*c(1,-1,-1,1)

n <- 100
eps <- 0.1
counts <- rmultinom(1, n, prob = joint_probs0)
priv_counts <- dp_count(counts, eps, delta = 2)
tab <- matrix(counts, 2, 2)
chisq.test(tab)

priv_tab <- matrix(priv_counts, 2, 2)
chisq.test(priv_tab)

fima_chi2(priv_tab, n, eps)

# Logistic
eps <- 1
H <- 10^4
n <- 1000
beta0 <- 0.5
beta1 <- 2
beta2 <- -2
beta <- c(beta0, beta1, beta2)
p1 <- 0.5
p2 <- 0.5

set.seed(123)
X1 <- rbinom(n, 1, p1)
X2 <- rbinom(n, 1, p2)
mu <- beta0 + beta1 * X1 + beta2 * X2
probs <- expit(mu)

# Approximation
W1 <- (X1 == 0) * (X2 == 0)
W2 <- (X1 == 1) * (X2 == 0)
W3 <- (X1 == 0) * (X2 == 1)
W4 <- (X1 == 1) * (X2 == 1)

set.seed(123)
u <- runif(n)
Y <- as.numeric(probs > u)

pi1 <- mean(Y[W1 == 1])
pi2 <- mean(Y[W2 == 1])
pi3 <- mean(Y[W3 == 1])
pi4 <- mean(Y[W4 == 1])

n1 <- sum(W1)
n2 <- sum(W2)
n3 <- sum(W3)
n4 <- sum(W4)

dp_pi <- dp_prop(c(pi1, pi2, pi3), eps = eps / 3, n = n)

dist <- fima_logistic(dp_pi, n, eps = 1, delta = 2, H = 10^4, seed = 123)
