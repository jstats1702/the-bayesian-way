# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
setwd("~/Dropbox/UN/bayes_book_en")

# Posterior distribution -------------------------------------------------------

# Model hyperparameters
a <- 1
b <- 1

# Observed data
n <- 20
y <- 1

# Parameters of the posterior distribution
a_post <- a + y
b_post <- b + n - y

# Posterior mean
media_post <- a_post / (a_post + b_post)

# Percentile-based interval ----------------------------------------------------

alpha <- 0.05

ic_perc <- qbeta(
     c(alpha / 2, 1 - alpha / 2),
     shape1 = a_post,
     shape2 = b_post
)

round(ic_perc, 4)

# HPD interval -----------------------------------------------------------------

# Function to compute the HPD interval for a unimodal Beta
hpd_beta <- function(shape1, shape2, prob = 0.95) {
     width <- function(l) {
          u <- qbeta(
               pbeta(l, shape1 = shape1, shape2 = shape2) + prob,
               shape1 = shape1,
               shape2 = shape2
          )
          
          u - l
     }
     
     upper_l <- qbeta(1 - prob, shape1 = shape1, shape2 = shape2)
     
     opt <- optimize(width, interval = c(0, upper_l))
     
     l <- opt$minimum
     u <- qbeta(
          pbeta(l, shape1 = shape1, shape2 = shape2) + prob,
          shape1 = shape1,
          shape2 = shape2
     )
     
     c(l, u)
}

# HPD interval
ic_hpd <- hpd_beta(
     shape1 = a_post,
     shape2 = b_post,
     prob   = 1 - alpha
)

round(ic_hpd, 4)

# Visualization ----------------------------------------------------------------

pdf(file = "prevalencia_intervalos.pdf", width = 5, height = 5, pointsize = 15)

par(mar = c(3.25, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

curve(
     expr = dbeta(x, shape1 = a_post, shape2 = b_post),
     from = 0,
     to   = 0.4,
     n    = 10000,
     col  = "black",
     lwd  = 1,
     xlab = expression(theta),
     ylab = expression(paste("p(", theta, " | ", y, ")")),
     main = "",
     xaxt = "n"
)

# Main x-axis
axis(side = 1, at = seq(0, 1, by = 0.2))

# Density height at the interval limits
y_perc <- dbeta(ic_perc, shape1 = a_post, shape2 = b_post)
y_hpd  <- dbeta(ic_hpd, shape1 = a_post, shape2 = b_post)

# Limits of the percentile-based interval
segments(
     x0  = ic_perc,
     y0  = 0,
     x1  = ic_perc,
     y1  = y_perc,
     col = "royalblue",
     lwd = 4,
     lty = 1
)

# Limits of the HPD interval
segments(
     x0  = ic_hpd,
     y0  = 0.2,
     x1  = ic_hpd,
     y1  = y_hpd,
     col = "darkgreen",
     lwd = 4,
     lty = 1
)

# Percentile-based interval
segments(
     x0  = ic_perc[1],
     y0  = 0,
     x1  = ic_perc[2],
     y1  = 0,
     col = "royalblue",
     lwd = 4
)

# HPD interval
segments(
     x0  = ic_hpd[1],
     y0  = 0.2,
     x1  = ic_hpd[2],
     y1  = 0.2,
     col = "darkgreen",
     lwd = 4
)

legend(
     "topright",
     legend = c("Posterior", "95% Percentile CI", "95% HPD CI"),
     col    = c("black", "royalblue", "darkgreen"),
     lty    = 1,
     lwd    = c(2, 2, 2),
     bty    = "n"
)

dev.off()

# End --------------------------------------------------------------------------