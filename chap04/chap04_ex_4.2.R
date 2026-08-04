# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
setwd("~/Dropbox/UN/bayes_book_en")

# Simulation parameters
set.seed(123)

a <- 3
b <- 2
B <- 1000

# Monte Carlo sample
theta_mc_3 <- rgamma(B, shape = a, rate = b)

# Summaries --------------------------------------------------------------------

# Mean
media_theta <- mean(theta_mc_3)

round(media_theta, 3)

# Standard deviation
desv_theta <- sd(theta_mc_3)

round(desv_theta, 3)

# Monte Carlo standard error
error_estandar <- desv_theta / sqrt(length(theta_mc_3))

round(error_estandar, 3)

# Monte Carlo coefficient of variation
cv_mc <- error_estandar / abs(media_theta)

round(cv_mc, 3)

# Sample size required for a margin of error of 0.01
epsilon     <- 0.01
B_necesario <- (1.96 * desv_theta / epsilon)^2

round(B_necesario)

# End --------------------------------------------------------------------------