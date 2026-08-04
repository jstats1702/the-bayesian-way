# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
setwd("~/Dropbox/UN/bayes_book_en")

# Data processing --------------------------------------------------------------

# Data
df <- read.delim("victimas.txt", stringsAsFactors = FALSE)

# Database dimensions
dim(df)

# Available variables
names(df)

# Sex frequencies
table(df$sexo, useNA = "ifany")

# Proportion with no sex information
round(mean(df$sexo == "Sin Informacion", na.rm = TRUE), 4)

# Sex coding
df <- df[df$sexo != "Sin Informacion", ]

df$sexo <- ifelse(df$sexo == "Mujer", 1, ifelse(df$sexo == "Hombre", 0, NA))

df$sexo <- as.numeric(df$sexo)

# Sex in 2016
y <- df[df$agno == 2016, "sexo"]

# Sex frequencies in 2016
table(y)

# Sample size
(n <- length(y))

# Sufficient statistic
(s <- sum(y))

# Modeling ---------------------------------------------------------------------

# Hyperparameters
a <- 1
b <- 1

# Posterior parameters
(ap <- a + s)
(bp <- b + n - s)

# Visualization
pdf(file = "victimas_posterior.pdf", width = 5, height = 5, pointsize = 15)

par(mar = c(3.25, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

# Plot of the posterior Beta distribution
curve(
     expr = dbeta(x, shape1 = ap, shape2 = bp),
     from = 0,
     to   = 1,
     n    = 10000,
     col  = "red",
     lwd  = 2,
     xlab = expression(theta),
     ylab = expression(paste("p(", theta, " | ", y, ")")),
     main = ""
)

# Beta(1,1) prior distribution
curve(
     expr = dbeta(x, shape1 = a, shape2 = b),
     from = 0,
     to   = 1,
     n    = 10000,
     col  = "royalblue",
     lwd  = 2,
     add  = TRUE
)

# Legend
legend(
     "topleft",
     legend = c("Prior", "Posterior"),
     col    = c("royalblue", "red"),
     lty    = 1,
     lwd    = 2,
     bty    = "n"
)

dev.off()

# Inference --------------------------------------------------------------------

# Computation of posterior summaries
media    <- ap / (ap + bp)
mediana  <- qbeta(0.5, shape1 = ap, shape2 = bp)
moda     <- (ap - 1) / (ap + bp - 2)
varianza <- (ap * bp) / ((ap + bp)^2 * (ap + bp + 1))
cv       <- sqrt(varianza) / media
ic_perc  <- qbeta(p = c(0.025, 0.975), shape1 = ap, shape2 = bp)

# Create results table
out <- data.frame(
     Media       = media,
     Mediana     = mediana,
     Moda        = moda,
     CV          = cv,
     `Q2.5%`     = ic_perc[1],
     `Q97.5%`    = ic_perc[2],
     check.names = FALSE
)

round(out, 4)

# Visualization
pdf(file = "victimas_intervalos.pdf", width = 5, height = 5, pointsize = 15)

par(mar = c(3.25, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

curve(
     expr = dbeta(x, shape1 = ap, shape2 = bp),
     from = 0.7,
     to   = 1,
     n    = 10000,
     col  = "red",
     lwd  = 2,
     xlab = expression(theta),
     ylab = expression(paste("p(", theta, " | ", y, ")")),
     main = ""
)

# 95% percentile-based credible interval
ic_perc <- qbeta(
     p      = c(0.025, 0.975),
     shape1 = ap,
     shape2 = bp
)

abline(
     v   = ic_perc,
     col = "darkgreen",
     lty = 4,
     lwd = 2
)

# Posterior mean
media <- ap / (ap + bp)

abline(
     v   = media,
     col = "black",
     lty = 3,
     lwd = 2
)

dev.off()

# End --------------------------------------------------------------------------