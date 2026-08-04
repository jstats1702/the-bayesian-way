# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
setwd("~/Dropbox/UN/bayes_book_en")

# Data ------------------------------------------------------------------------

n      <- 1000
y      <- 532
theta0 <- 0.5

# Frequentist approach ---------------------------------------------------------

# Test statistic using the Normal approximation
z <- (y - n * theta0) / sqrt(n * theta0 * (1 - theta0))

round(z, 3)

# Approximate two-sided p-value
p_val <- 2 * pnorm(q = abs(z), lower.tail = FALSE)

round(p_val, 3)

# Bayesian approach ------------------------------------------------------------

# Prior probabilities of the hypotheses
p_H0 <- 0.5
p_H1 <- 0.5

# Under H0: theta = 0.5
py_H0 <- dbinom(x = y, size = n, prob = theta0)

round(py_H0, 4)

# Under H1: theta | H1 ~ Beta(a, b)
a <- 1
b <- 1

# Marginal likelihood under H1 using the Beta-Binomial distribution
py_H1 <- exp(
     lchoose(n, y) +
          lbeta(a + y, b + n - y) -
          lbeta(a, b)
)

round(py_H1, 4)

# Equivalent expression for a = b = 1
round(1 / (n + 1), 4)

# Bayes factor in favor of H0 over H1
B01 <- py_H0 / py_H1

round(B01, 3)

# Bayes factor in favor of H1 over H0
B10 <- 1 / B01

round(B10, 3)

# Posterior probabilities of the hypotheses -----------------------------------

post_H0 <- (py_H0 * p_H0) / (py_H0 * p_H0 + py_H1 * p_H1)
post_H1 <- (py_H1 * p_H1) / (py_H0 * p_H0 + py_H1 * p_H1)

round(post_H0, 3)
round(post_H1, 3)

# End --------------------------------------------------------------------------