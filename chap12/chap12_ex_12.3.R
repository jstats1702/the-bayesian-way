# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
setwd("~/Dropbox/UN/bayes_book_en")

# Load the required package
suppressMessages(suppressWarnings(library(stats)))

# Parameters
alpha <- 10
eps   <- 1e-6

# Select the truncation level
H <- ceiling(
     log(eps) /
          log(alpha / (alpha + 1))
)

# Base distribution
G0 <- function(x) {
     pnorm(x)
}

set.seed(123)

# Simulate the atoms
vartheta <- rnorm(
     n    = H,
     mean = 0,
     sd   = 1
)

# Simulate the sequential stick-breaking process variables
V <- rbeta(
     n      = H,
     shape1 = 1,
     shape2 = alpha
)

# Compute the weights before reallocating the residual mass
omega <- V * c(
     1,
     cumprod(1 - V[-H])
)

# Compute the residual mass after H breaks
R_H <- prod(1 - V)

# Assign the residual mass to the last weight
omega[H] <- omega[H] + R_H

# Evaluate the distribution function of the realization
x_vals <- seq(
     from       = -3,
     to         = 3,
     length.out = 10000
)

G <- sapply(
     X   = x_vals,
     FUN = function(x) {
          sum(omega[vartheta <= x])
     }
)

# Visualization of the weights
pdf(
     file      = "dp_particion_secuencial_pesos.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Left panel: atoms and weights
plot(
     x    = vartheta,
     y    = omega,
     type = "h",
     xlim = c(-3, 3),
     xlab = expression(vartheta),
     ylab = expression(omega),
     main = "",
     col  = 4,
     lend = 1
)

dev.off()

# Visualization of the realization
pdf(
     file      = "dp_particion_secuencial_realizacion.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

# Configure the plotting window
par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Right panel: distribution function of the realization
plot(
     x    = x_vals,
     y    = G,
     type = "l",
     xlim = c(-3, 3),
     ylim = c(0, 1),
     xlab = "x",
     ylab = "G(x)",
     main = "",
     col  = 4
)

# Add the base distribution function
curve(
     expr = G0,
     from = -3,
     to   = 3,
     n    = 1000,
     lwd  = 2,
     add  = TRUE
)

dev.off()

# End --------------------------------------------------------------------------