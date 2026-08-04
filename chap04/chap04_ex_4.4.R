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

# Observed statistics ----------------------------------------------------------

t_obs_1 <- c(media = mean(y1), de = sd(y1))
t_obs_2 <- c(media = mean(y2), de = sd(y2))

round(t_obs_1, 3)
round(t_obs_2, 3)

# Monte Carlo simulation -------------------------------------------------------

# Number of samples
B <- 10000

# Posterior samples of theta_1 and theta_2
set.seed(123)
th1_mc <- rgamma(B, shape = ap1, rate = bp1)
th2_mc <- rgamma(B, shape = ap2, rate = bp2)

# Initialize matrices to store test statistics
set.seed(123)
t_mc_1 <- matrix(NA, nrow = B, ncol = 2)
t_mc_2 <- matrix(NA, nrow = B, ncol = 2)

colnames(t_mc_1) <- c("media", "de")
colnames(t_mc_2) <- c("media", "de")

# Posterior predictive distribution
set.seed(123)

for (i in seq_len(B)) {
     # Replicated data
     y1_rep <- rpois(n = n1, lambda = th1_mc[i])
     y2_rep <- rpois(n = n2, lambda = th2_mc[i])
     
     # Test statistics
     t_mc_1[i, ] <- c(mean(y1_rep), sd(y1_rep))
     t_mc_2[i, ] <- c(mean(y2_rep), sd(y2_rep))
}

# Posterior predictive p-values -----------------------------------------------

ppp_1 <- c(
     media = mean(t_mc_1[, "media"] <= t_obs_1["media"]),
     de    = mean(t_mc_1[, "de"] <= t_obs_1["de"])
)

ppp_2 <- c(
     media = mean(t_mc_2[, "media"] <= t_obs_2["media"]),
     de    = mean(t_mc_2[, "de"] <= t_obs_2["de"])
)

round(ppp_1, 3)
round(ppp_2, 3)

# Visualization ----------------------------------------------------------------

# Colors
col1 <- 2
col2 <- 4

# Axis limits
xlim_media <- range(
     t_mc_1[, "media"],
     t_mc_2[, "media"],
     t_obs_1["media"],
     t_obs_2["media"]
)

xlim_de <- range(
     t_mc_1[, "de"],
     t_mc_2[, "de"],
     t_obs_1["de"],
     t_obs_2["de"]
)

ylim_media <- c(0, 4)
ylim_de    <- c(0, 6.5)

# Histogram of the mean, group without higher education ------------------------

pdf(file = "hijos_chequeo_media_sin.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

hist(
     x      = t_mc_1[, "media"],
     freq   = FALSE,
     col    = col1,
     border = col1,
     xlim   = xlim_media,
     ylim   = ylim_media,
     xlab   = "t",
     ylab   = expression(p(t ~ "|" ~ bold(y))),
     main   = ""
)

abline(v = t_obs_1["media"], col = 1, lwd = 2)

legend(
     "topleft",
     legend = c("Without", "With", "Obs"),
     bty    = "n",
     fill   = c(2, 4, 1),
     border = c(2, 4, 1),
     col    = c(2, 4, 1)
)

dev.off()

# Histogram of the mean, group with higher education ---------------------------

pdf(file = "hijos_chequeo_media_con.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

hist(
     x      = t_mc_2[, "media"],
     freq   = FALSE,
     col    = col2,
     border = col2,
     xlim   = xlim_media,
     ylim   = ylim_media,
     xlab   = "t",
     ylab   = expression(p(t ~ "|" ~ bold(y))),
     main   = ""
)

abline(v = t_obs_2["media"], col = 1, lwd = 2)

dev.off()

# Histogram of the standard deviation, group without higher education ----------

pdf(file = "hijos_chequeo_de_sin.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

hist(
     x      = t_mc_1[, "de"],
     freq   = FALSE,
     col    = col1,
     border = col1,
     xlim   = xlim_de,
     ylim   = ylim_de,
     xlab   = "t",
     ylab   = expression(p(t ~ "|" ~ bold(y))),
     main   = ""
)

abline(v = t_obs_1["de"], col = 1, lwd = 2)

legend(
     "topleft",
     legend = c("Without", "With", "Obs"),
     bty    = "n",
     fill   = c(2, 4, 1),
     border = c(2, 4, 1),
     col    = c(2, 4, 1)
)

dev.off()

# Histogram of the standard deviation, group with higher education -------------

pdf(file = "hijos_chequeo_de_con.pdf", width = 5, height = 5, pointsize = 15)

par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))

hist(
     x      = t_mc_2[, "de"],
     freq   = FALSE,
     col    = col2,
     border = col2,
     xlim   = xlim_de,
     ylim   = ylim_de,
     xlab   = "t",
     ylab   = expression(p(t ~ "|" ~ bold(y))),
     main   = ""
)

abline(v = t_obs_2["de"], col = 1, lwd = 2)

dev.off()

# End --------------------------------------------------------------------------