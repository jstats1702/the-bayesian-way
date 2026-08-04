# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
setwd("~/Dropbox/UN/bayes_book_en")

# Libraries
suppressMessages(suppressWarnings(library(stats)))

alpha   <- 5
n_small <- 10
n_large <- 100

x_grid <- seq(
     from       = -5,
     to         = 5,
     length.out = 2000
)

# Sample generation ------------------------------------------------------------

set.seed(123)

y_large <- rcauchy(
     n        = n_large,
     location = 0,
     scale    = 1
)

# The small sample corresponds to the first ten observations
y_small <- y_large[seq_len(n_small)]

# Distribution functions -------------------------------------------------------

# True distribution function
F_true <- function(x) {
     pcauchy(
          q        = x,
          location = 0,
          scale    = 1
     )
}

# Prior base distribution
G0 <- function(x) {
     pnorm(
          q    = x,
          mean = 0,
          sd   = 1
     )
}

# Empirical distribution functions
G_emp_small <- ecdf(y_small)
G_emp_large <- ecdf(y_large)

# Posterior means
G_post_small <- function(x) {
     (alpha * G0(x) + n_small * G_emp_small(x)) / (alpha + n_small)
}

G_post_large <- function(x) {
     (alpha * G0(x) + n_large * G_emp_large(x)) / (alpha + n_large)
}

# Weights of the prior and empirical distributions -----------------------------

weights <- data.frame(
     n                = c(n_small, n_large),
     weight_G0        = alpha / (alpha + c(n_small, n_large)),
     weight_empirical = c(n_small, n_large) /
          (alpha + c(n_small, n_large))
)

round(weights, 3)

# Function for constructing each plot ------------------------------------------

plot_cdf_comparison <- function(
          G_emp,
          G_post,
          show_legend = FALSE
) {
     plot(
          x    = x_grid,
          y    = F_true(x_grid),
          type = "l",
          lty  = 1,
          lwd  = 2.5,
          col  = "black",
          ylim = c(0, 1),
          xlab = expression(x),
          ylab = "Distribution function",
          main = ""
     )
     
     lines(
          x    = x_grid,
          y    = G0(x_grid),
          type = "l",
          lty  = 1,
          lwd  = 2,
          col  = "gray60"
     )
     
     lines(
          x    = x_grid,
          y    = G_emp(x_grid),
          type = "s",
          lty  = 1,
          lwd  = 1.75,
          col  = "royalblue3"
     )
     
     lines(
          x    = x_grid,
          y    = G_post(x_grid),
          type = "l",
          lty  = 1,
          lwd  = 2.5,
          col  = "firebrick3"
     )
     
     if (show_legend) {
          legend(
               "topleft",
               legend = c(
                    "Cauchy(0, 1)",
                    "G0 = N(0, 1)",
                    "Empirical distribution",
                    "Posterior mean"
               ),
               col = c(
                    "black",
                    "gray60",
                    "royalblue3",
                    "firebrick3"
               ),
               lty = 1,
               lwd = c(2.5, 2, 1.75, 2.5),
               bty = "n",
               cex = 0.85
          )
     }
}

# Figure for n = 10 ------------------------------------------------------------

pdf(
     file      = "dp_cauchy_cdf_posterior_n10.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot_cdf_comparison(
     G_emp       = G_emp_small,
     G_post      = G_post_small,
     show_legend = TRUE
)

dev.off()

# Figure for n = 100 -----------------------------------------------------------

pdf(
     file      = "dp_cauchy_cdf_posterior_n100.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot_cdf_comparison(
     G_emp       = G_emp_large,
     G_post      = G_post_large,
     show_legend = FALSE
)

dev.off()

# End --------------------------------------------------------------------------