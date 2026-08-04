# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
setwd("~/Dropbox/UN/bayes_book_en")

# Data ------------------------------------------------------------------------

# Data: length of stay
LoS <- c(1, 2, 1, 1, 4, 1, 2, 2, 0, 3, 6, 2, 1, 3)

(n <- length(LoS))
(y <- sum(LoS))

# Prior distribution ----------------------------------------------------------

a <- 1
b <- 1 / 3

# Prior mean
round(a / b, 3)

# Prior coefficient of variation
round(1 / sqrt(a), 3)

# Posterior distribution -------------------------------------------------------

# Posterior parameters
(ap <- a + y)
(bp <- b + n)

# Inference -------------------------------------------------------------------

# Posterior mean
round(ap / bp, 3)

# Posterior coefficient of variation
round(1 / sqrt(ap), 3)

# 95% credible interval
round(qgamma(p = c(0.025, 0.975), shape = ap, rate = bp), 3)

# One-sided test: H0: theta >= 3 vs H1: theta < 3 -----------------------------

theta0 <- 3

# Posterior probabilities of the hypotheses
post_H0 <- pgamma(q = theta0, shape = ap, rate = bp, lower.tail = FALSE)
post_H1 <- pgamma(q = theta0, shape = ap, rate = bp, lower.tail = TRUE)

# P(H0 | y)
round(post_H0, 3)

# P(H1 | y)
round(post_H1, 3)

# Prior probabilities of the hypotheses
prior_H0 <- pgamma(q = theta0, shape = a, rate = b, lower.tail = FALSE)
prior_H1 <- pgamma(q = theta0, shape = a, rate = b, lower.tail = TRUE)

# P(H0)
round(prior_H0, 3)

# P(H1)
round(prior_H1, 3)

# Bayes factor in favor of H1 over H0
B10 <- (post_H1 / post_H0) / (prior_H1 / prior_H0)

round(B10, 3)

# Bayes factor in favor of H0 over H1
B01 <- 1 / B10

round(B01, 3)

# End --------------------------------------------------------------------------