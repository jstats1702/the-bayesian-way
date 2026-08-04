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
     df$P_NIVEL_ANOSR %in% c(8, 9),
     1,
     ifelse(is.na(df$P_NIVEL_ANOSR) | df$P_NIVEL_ANOSR == 99, NA, 0)
)

# Frequencies: higher education indicator
table(df$educ_sup, useNA = "ifany")

round(100 * prop.table(table(df$educ_sup, useNA = "ifany")), 1)

# Live-born children -----------------------------------------------------------

# PA1_THNV: Total number of live-born children
table(df$PA1_THNV, useNA = "ifany")

# Recoding: NA is interpreted as 0 children due to questionnaire skip logic
df$hijos <- as.numeric(replace(df$PA1_THNV, is.na(df$PA1_THNV), 0))

# Frequencies: number of children
table(df$hijos, useNA = "ifany")

# 0 children
sum(df$hijos == 0)
round(100 * sum(df$hijos == 0) / length(df$hijos), 3)

# 1 or more children
sum((df$hijos != 0) & (df$hijos != 99))
round(100 * sum((df$hijos != 0) & (df$hijos != 99)) / length(df$hijos), 3)

# Number of children not reported
sum(df$hijos == 99)
round(100 * sum(df$hijos == 99) / length(df$hijos), 3)

# Remove missing or unreported data -------------------------------------------

df <- subset(df, !is.na(educ_sup) & hijos != 99)

# Selection filter -------------------------------------------------------------

filtro <- with(
     df,
     (P_PARENTESCOR == 1) &
          (P_SEXO == 2) &
          (P_EDADR == 9) &
          (PA1_GRP_ETNIC == 6) &
          (PA_LUG_NAC %in% c(2, 3)) &
          (PA_VIVIA_5ANOS %in% c(2, 3)) &
          (PA_HNV %in% c(1, 2)) &
          (P_ALFABETA == 1)
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

# Frequency distribution -------------------------------------------------------

pdf(file = "hijos_barras.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

y_vals <- 0:6

freq_y1 <- table(factor(y1, levels = y_vals)) / n1
freq_y2 <- table(factor(y2, levels = y_vals)) / n2

freq_mat <- rbind(
     "Sin" = freq_y1,
     "Con" = freq_y2
)

barplot(
     freq_mat,
     beside      = TRUE,
     ylim        = c(0, 0.4),
     names.arg   = y_vals,
     ylab        = "Relative frequency",
     xlab        = "No. of children",
     legend.text = rownames(freq_mat),
     args.legend = list(
          bty    = "n",
          x      = "topright",
          border = c(2, 4)
     ),
     col         = c(2, 4),
     border      = NA
)

dev.off()

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

# Visualization ----------------------------------------------------------------

pdf(file = "hijos_posterior.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

# Prior distribution and posterior distributions
theta <- seq(0, 5, length.out = 1000)

plot(
     NA,
     NA,
     xlim = c(0, 4),
     ylim = c(0, 5.5),
     xlab = expression(theta),
     ylab = expression(paste("p(", theta, " | ", y, ")")),
     main = ""
)

lines(theta, dgamma(theta, shape = ap1, rate = bp1), col = 2, lwd = 2)
lines(theta, dgamma(theta, shape = ap2, rate = bp2), col = 4, lwd = 2)
lines(theta, dgamma(theta, shape = a, rate = b), col = 1, lwd = 1)

abline(h = 0, col = 1)

legend(
     "topright",
     legend = c("Without", "With", "Prior"),
     bty    = "n",
     lwd    = c(2, 2, 1),
     col    = c(2, 4, 1)
)

dev.off()

# Inference -------------------------------------------------------------------

# Posterior mean of theta
theta_hat_1 <- ap1 / bp1
theta_hat_2 <- ap2 / bp2

# Posterior coefficient of variation of theta
cv_1 <- 1 / sqrt(ap1)
cv_2 <- 1 / sqrt(ap2)

# 95% credible interval for theta
ic95_1 <- qgamma(p = c(0.025, 0.975), shape = ap1, rate = bp1)
ic95_2 <- qgamma(p = c(0.025, 0.975), shape = ap2, rate = bp2)

# Posterior probability that theta > 2
pr_theta_1 <- pgamma(q = 2, shape = ap1, rate = bp1, lower.tail = FALSE)
pr_theta_2 <- pgamma(q = 2, shape = ap2, rate = bp2, lower.tail = FALSE)

# Results table ---------------------------------------------------------------

tab <- rbind(
     "Sin superior" = c(theta_hat_1, cv_1, ic95_1, pr_theta_1),
     "Con superior" = c(theta_hat_2, cv_2, ic95_2, pr_theta_2)
)

colnames(tab) <- c("Media", "CV", "Q2.5%", "Q97.5%", "Pr. > 2")

round(tab, 3)

# Hypothesis testing -----------------------------------------------------------

# H0: theta1 = theta2, with theta ~ Gamma(a1, b1)
# H1: theta1 ~ Gamma(a1, b1), theta2 ~ Gamma(a2, b2)

# Hyperparameters under H1: theta1 ~ Gamma(a1, b1), theta2 ~ Gamma(a2, b2)
a1 <- a
b1 <- b
a2 <- a
b2 <- b

# Bayes factor B01 on the logarithmic scale
log_B01 <- (
     a1 * log(b1) - lgamma(a1) +
          lgamma(a1 + s1 + s2) -
          (a1 + s1 + s2) * log(b1 + n1 + n2)
) - (
     a1 * log(b1) - lgamma(a1) +
          lgamma(a1 + s1) -
          (a1 + s1) * log(b1 + n1) +
          a2 * log(b2) - lgamma(a2) +
          lgamma(a2 + s2) -
          (a2 + s2) * log(b2 + n2)
)

# Bayes factors
B01 <- exp(log_B01)
B10 <- exp(-log_B01)

round(B01, 6)
round(B10, 3)

# Posterior probabilities, assuming P(H0) = P(H1) = 0.5
post_H1 <- B10 / (1 + B10)
post_H0 <- 1 / (1 + B10)

round(post_H1, 6)
round(post_H0, 6)

# Frequentist p-value ----------------------------------------------------------

# p-value = Pr(observing a difference as extreme or more extreme under H0)
yb1 <- mean(y1)
yb2 <- mean(y2)

sd1 <- sd(y1)
sd2 <- sd(y2)

z_welch <- (yb1 - yb2) / sqrt(sd1^2 / n1 + sd2^2 / n2)
p_welch <- 2 * pnorm(q = abs(z_welch), lower.tail = FALSE)

round(z_welch, 3)
p_welch

# End --------------------------------------------------------------------------