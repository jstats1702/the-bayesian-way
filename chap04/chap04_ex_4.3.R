# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
setwd("~/Dropbox/UN/bayes_book_en")

# Data ------------------------------------------------------------------------

df <- read.csv("CNPV2018.txt", stringsAsFactors = FALSE)

# Dimensions
dim(df)

# Recoding of educational attainment ------------------------------------------

# 0: No higher education
# 1: Higher education
df$educ_sup <- ifelse(
     df$P_NIVEL_ANOSR %in% c(8, 9), 1,
     ifelse(is.na(df$P_NIVEL_ANOSR) | df$P_NIVEL_ANOSR == 99, NA, 0)
)

# Frequencies of the higher education indicator
table(df$educ_sup, useNA = "ifany")
round(100 * prop.table(table(df$educ_sup, useNA = "ifany")), 1)

# Live-born children -----------------------------------------------------------

# PA1_THNV: total number of live-born children
table(df$PA1_THNV, useNA = "ifany")

# Recoding: NA is interpreted as 0 children due to questionnaire skip logic
df$hijos <- as.numeric(replace(df$PA1_THNV, is.na(df$PA1_THNV), 0))

# Frequencies of the number of children
table(df$hijos, useNA = "ifany")

# 0 children
sum(df$hijos == 0)
round(100 * sum(df$hijos == 0) / length(df$hijos), 3)

# 1 or more children
sum(df$hijos != 0 & df$hijos != 99)
round(100 * sum(df$hijos != 0 & df$hijos != 99) / length(df$hijos), 3)

# Number of children not reported
sum(df$hijos == 99)
round(100 * sum(df$hijos == 99) / length(df$hijos), 3)

# Remove missing or unreported data -------------------------------------------

df <- subset(df, !is.na(educ_sup) & hijos != 99)

# Selection filter -------------------------------------------------------------

filtro <- with(
     df,
     P_PARENTESCOR == 1 &
          P_SEXO == 2 &
          P_EDADR == 9 &
          PA1_GRP_ETNIC == 6 &
          PA_LUG_NAC %in% c(2, 3) &
          PA_VIVIA_5ANOS %in% c(2, 3) &
          PA_HNV %in% c(1, 2) &
          P_ALFABETA == 1
)

# Number of children by educational attainment --------------------------------

y1 <- df$hijos[filtro & df$educ_sup == 0]  # No higher education
y2 <- df$hijos[filtro & df$educ_sup == 1]  # Higher education

# Sample sizes
(n1 <- length(y1))
(n2 <- length(y2))

# Sufficient statistics
(s1 <- sum(y1))
(s2 <- sum(y2))

# Gamma(2, 1) prior ------------------------------------------------------------

a <- 2
b <- 1

# Prior mean of theta
round(a / b, 3)

# Prior coefficient of variation of theta
round(1 / sqrt(a), 3)

# Posterior parameters ---------------------------------------------------------

# No higher education
(ap1 <- a + s1)
(bp1 <- b + n1)

# Higher education
(ap2 <- a + s2)
(bp2 <- b + n2)

# Monte Carlo simulation -------------------------------------------------------

# Number of samples
B <- 10000

# Posterior samples of theta_1 and theta_2
set.seed(123)
th1_mc <- rgamma(B, shape = ap1, rate = bp1)
th2_mc <- rgamma(B, shape = ap2, rate = bp2)

# Posterior sample of Delta = theta_1 - theta_2
Delta <- th1_mc - th2_mc

# Posterior summaries of Delta -----------------------------------------------

# Posterior mean of Delta
round(mean(Delta), 3)

# Posterior coefficient of variation of Delta
round(sd(Delta) / abs(mean(Delta)), 3)

# 95% credible interval for Delta
round(quantile(Delta, probs = c(0.025, 0.975)), 3)

# Posterior probability that Delta > 0
round(mean(Delta > 0), 3)

# Plots ------------------------------------------------------------------------

pdf(file = "hijos_diferencia.pdf", width = 5, height = 5, pointsize = 15)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

hist(
     Delta,
     breaks = 30,
     freq   = FALSE,
     col    = "gold2",
     border = "gold2",
     xlab   = expression(Delta),
     ylab   = expression(paste("p(", Delta, " | ", y, ")")),
     main   = ""
)

lines(density(Delta), col = "black", lwd = 2)

dev.off()

# Posterior predictive simulation ---------------------------------------------

# Samples from the posterior predictive distribution
set.seed(123)
y1_mc <- rpois(B, lambda = th1_mc)
y2_mc <- rpois(B, lambda = th2_mc)

# Posterior predictive sample of d = y*_1 - y*_2
d_mc <- y1_mc - y2_mc

# Posterior predictive mean of d
round(mean(d_mc), 3)

# 95% posterior predictive interval for d
round(quantile(d_mc, probs = c(0.025, 0.975)), 3)

# Posterior predictive probability that d > 0
round(mean(d_mc > 0), 3)

# Plot of the posterior predictive distributions ------------------------------

pdf(file = "hijos_predictiva.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(2, 0.75, 0))

y_vals <- 0:6

pred_y1 <- dnbinom(y_vals, size = ap1, mu = ap1 / bp1)
pred_y2 <- dnbinom(y_vals, size = ap2, mu = ap2 / bp2)

pred_mat <- rbind(
     "Sin" = pred_y1,
     "Con" = pred_y2
)

barplot(
     pred_mat,
     beside      = TRUE,
     ylim        = c(0, 0.4),
     names.arg   = y_vals,
     xlab        = expression(y^"*"),
     ylab        = expression(p(y^"*" ~ "|" ~ bold(y))),
     main        = "",
     legend.text = rownames(pred_mat),
     args.legend = list(
          bty    = "n",
          x      = "topright",
          border = c(2, 4)
     ),
     col         = c(2, 4),
     border      = NA
)

dev.off()

# Plot of the posterior predictive difference ---------------------------------

pdf(file = "hijos_predictiva_diferencia.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(2, 0.75, 0))

d_vals <- min(d_mc):max(d_mc)

d_freq <- table(factor(d_mc, levels = d_vals)) / B

barplot(
     d_freq,
     ylim      = c(0, max(d_freq) * 1.1),
     names.arg = d_vals,
     xlab      = expression(d == y[1]^"*" - y[2]^"*"),
     ylab      = expression(p(d ~ "|" ~ bold(y))),
     main      = "",
     col       = "darkgrey",
     border    = NA
)

dev.off()

# End --------------------------------------------------------------------------