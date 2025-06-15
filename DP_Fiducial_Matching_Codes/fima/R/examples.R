
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
