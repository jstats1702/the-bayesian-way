# Settings ---------------------------------------------------------------------

rm(list = ls())

setwd("~/Dropbox/UN/bayes_book_en")

# Data -------------------------------------------------------------------------

# Simulate data
n      <- 5
sigma2 <- 1

set.seed(123)

y <- rnorm(
     n    = n,
     mean = 10,
     sd   = sqrt(sigma2)
)

# Prior
mu0  <- 5
tau2 <- 10

# Log-posterior kernel ----------------------------------------------------------

log_posterior <- function(theta, y, sigma2, mu0, tau2) {
     log_likelihood <- sum(
          dnorm(
               x    = y,
               mean = theta,
               sd   = sqrt(sigma2),
               log  = TRUE
          )
     )
     
     log_prior <- dnorm(
          x    = theta,
          mean = mu0,
          sd   = sqrt(tau2),
          log  = TRUE
     )
     
     log_likelihood + log_prior
}

# Metropolis algorithm ---------------------------------------------------------

B      <- 10000  # number of iterations
theta  <- 0      # initial value
delta2 <- 2      # proposal variance
acept  <- 0      # acceptance counter

THETA <- numeric(B)  # storage

set.seed(123)

for (b in seq_len(B)) {
     # 1. Generate proposal
     theta_star <- rnorm(
          n    = 1,
          mean = theta,
          sd   = sqrt(delta2)
     )
     
     # 2. Compute the acceptance ratio on the logarithmic scale
     log_post_star <- log_posterior(theta_star, y, sigma2, mu0, tau2)
     log_post_curr <- log_posterior(theta, y, sigma2, mu0, tau2)
     
     log_r <- log_post_star - log_post_curr
     
     # 3. Accept or reject the proposal
     if (log(runif(1)) <= min(0, log_r)) {
          theta <- theta_star
          acept <- acept + 1
     }
     
     # 4. Store current state
     THETA[b] <- theta
}

# Acceptance rate --------------------------------------------------------------

tasa_aceptacion <- acept / B
round(tasa_aceptacion, 3)

# Chain ------------------------------------------------------------------------

# Trace plot
pdf(
     file      = "metropolis_ejemplo_juguete_traceplot.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = seq_along(THETA),
     y    = THETA,
     type = "p",
     pch  = 16,
     cex  = 0.25,
     col  = adjustcolor("black", 0.5),
     xlab = "Iteration",
     ylab = expression(theta)
)

dev.off()

# Posterior histogram ----------------------------------------------------------

pdf(
     file      = "metropolis_ejemplo_juguete_histograma.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Posterior samples without the warm-up period
theta_post <- THETA[-seq_len(50)]

hist(
     x      = theta_post,
     freq   = FALSE,
     col    = "gray80",
     border = "white",
     main   = "",
     xlab   = expression(theta),
     ylab   = "Density"
)

# Analytical posterior
theta_grid <- seq(
     from       = min(theta_post),
     to         = max(theta_post),
     length.out = 1000
)

mu_n <- (mu0 / tau2 + n * mean(y) / sigma2) /
     (1 / tau2 + n / sigma2)

tau2_n <- 1 / (1 / tau2 + n / sigma2)

lines(
     x   = theta_grid,
     y   = dnorm(theta_grid, mean = mu_n, sd = sqrt(tau2_n)),
     col = 2,
     lty = 1,
     lwd = 2
)

box()

dev.off()

# Sensitivity to the tuning parameter ------------------------------------------

# Proposal variance values
delta2_grid <- 2^c(-5, 1, 7)

# Algorithm configuration
B       <- 10000
burn_in <- 10

# Storage
n_delta <- length(delta2_grid)

ACR <- numeric(n_delta)  # acceptance rates
ACF <- numeric(n_delta)  # lag-1 autocorrelations
SAM <- matrix(NA, nrow = B, ncol = n_delta)

colnames(SAM) <- paste0("delta2_", delta2_grid)

for (k in seq_along(delta2_grid)) {
     delta2 <- delta2_grid[k]
     
     # Initial state
     theta <- 0
     acept <- 0
     THETA <- numeric(B)
     
     set.seed(123)
     
     for (b in seq_len(B)) {
          # 1. Generate proposal
          theta_star <- rnorm(
               n    = 1,
               mean = theta,
               sd   = sqrt(delta2)
          )
          
          # 2. Compute the acceptance ratio on the logarithmic scale
          log_post_star <- log_posterior(theta_star, y, sigma2, mu0, tau2)
          log_post_curr <- log_posterior(theta, y, sigma2, mu0, tau2)
          
          log_r <- log_post_star - log_post_curr
          
          # 3. Accept or reject the proposal
          if (log(runif(1)) <= min(0, log_r)) {
               theta <- theta_star
               acept <- acept + 1
          }
          
          # 4. Store current state
          THETA[b] <- theta
     }
     
     # Store results
     ACR[k]   <- acept / B
     ACF[k]   <- acf(THETA[-seq_len(burn_in)], plot = FALSE)$acf[2]
     SAM[, k] <- THETA
}

# Summary ----------------------------------------------------------------------

resumen_sensibilidad <- data.frame(
     delta2         = round(delta2_grid, 3),
     tasa_aceptacion = round(ACR, 3),
     acf_orden_1     = round(ACF, 3)
)

resumen_sensibilidad

# Plots ------------------------------------------------------------------------

for (k in seq_along(delta2_grid)) {
     pdf(
          file      = paste0(
               "metropolis_ejemplo_juguete_traceplot_delta_",
               k,
               ".pdf"
          ),
          width     = 5,
          height    = 5,
          pointsize = 17
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     plot(
          x    = seq_len(500),
          y    = SAM[seq_len(500), k],
          type = "p",
          pch  = 16,
          cex  = 0.4,
          col  = adjustcolor("black", 0.5),
          xlab = "Iteration",
          ylab = expression(theta),
          ylim = range(SAM)
     )
     
     abline(
          h   = mu_n,
          col = 2,
          lwd = 2
     )
     
     dev.off()
}

# End --------------------------------------------------------------------------