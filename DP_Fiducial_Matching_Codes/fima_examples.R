# One prop confidence interval
p <- 0.8
eps <- 1
alpha <- 0.05
H <- 10^4

n <- 30
set.seed(14)
x <- rbinom(n, 1, prob = p)

pi <- dp_prop(mean(x), eps = eps, n = n)

dist <- fima_prop(pi, eps = eps, n = n, H = H)
fima_infer(dist)


# Two prop hypothesis test

p1 <- 0.8
p2 <- 0.9
n <- 300

set.seed(14)
x1 <- rbinom(n, 1, prob = p1)
x2 <- rbinom(n, 1, prob = p2)

pi1 <- dp_prop(mean(x1), eps = 1, n = n)
pi2 <- dp_prop(mean(x2), eps = 1, n = n)

dist <- fima_2prop(pi1, pi2, n1 = n, n2 = n, eps = eps, H = H)
fima_infer(dist, theta = 0, ha = "smaller")


# Chi-squared
diff <- 0.01
joint_probs0 <- c(0.25, 0.25, 0.25, 0.25) #+ diff*c(1,-1,-1,1)

n <- 5000
eps <- 1

set.seed(124)
counts <- rmultinom(1, n, prob = joint_probs0)
tab <- matrix(counts, 2, 2)
chisq.test(tab)

priv_counts <- dp_count(counts, eps, delta = 2)
priv_tab <- matrix(priv_counts, 2, 2)
dist <- fima_chi2(priv_tab, n, eps)
fima_infer(dist)

# Logistic (this needs to be fixed!)
eps <- 1
H <- 10^4
n <- 1000
beta0 <- 0.5
beta1 <- 2
beta2 <- -2
p1 <- 0.5
p2 <- 0.5

set.seed(12)

X1 <- rbinom(n, 1, p1)
X2 <- rbinom(n, 1, p2)
mu <- beta0 + beta1 * X1 + beta2 * X2
probs <- expit(mu)

# Approximation
W1 <- (X1 == 0) * (X2 == 0)
W2 <- (X1 == 1) * (X2 == 0)
W3 <- (X1 == 0) * (X2 == 1)
W4 <- (X1 == 1) * (X2 == 1)

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

dp_pi <- c(dp_prop(pi1, eps = eps / 4, n = n1),
           dp_prop(pi2, eps = eps / 4, n = n2),
           dp_prop(pi3, eps = eps / 4, n = n3),
           dp_prop(pi4, eps = eps / 4, n = n4))

dist <- fima_logistic(dp_pi, n, eps = 1, delta = 2, H = 10^4, seed = 1234)
fima_infer(dist)
