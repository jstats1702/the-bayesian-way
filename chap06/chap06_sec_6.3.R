# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
setwd("~/Dropbox/UN/bayes_book_en")

# Load data
dat <- read.csv(
     file             = file.path("Examen_Saber_11_20251.txt"),
     sep              = ";",
     stringsAsFactors = FALSE
)

dim(dat)

# Initial review by department
table(dat$estu_depto_reside, useNA = "ifany")

round(100 * table(dat$estu_depto_reside, useNA = "ifany") / nrow(dat), 1)

# Clean department variable
dat$estu_depto_reside <- trimws(dat$estu_depto_reside)

# Remove records without a department and records from abroad
dat <- dat[
     !is.na(dat$estu_depto_reside) &
          dat$estu_depto_reside != "" &
          !dat$estu_depto_reside %in% c("EXTRANJERO", "EXTRANGERO"),
     ,
     drop = FALSE
]

# Verify remaining departments
table(dat$estu_depto_reside)
unique(dat$estu_depto_reside)

# Mathematics scores for Bogotá
y_bogota <- dat[
     dat$estu_depto_reside == "BOGOTÁ",
     "punt_matematicas"
]

# Remove missing values
y_bogota <- as.numeric(y_bogota)
y_bogota <- y_bogota[is.finite(y_bogota)]

# Sample size and mean
length(y_bogota)
round(mean(y_bogota), 3)

# Proportional sampling by department -----------------------------------------

prop_muestra  <- 0.05
departamentos <- sort(unique(dat$estu_depto_reside))

set.seed(123)

idx_muestra <- integer(0)

for (depto in departamentos) {
     idx_depto     <- which(dat$estu_depto_reside == depto)
     n_depto_total <- length(idx_depto)
     
     if (n_depto_total < 100) {
          n_depto <- min(10, n_depto_total)
     } else {
          n_depto <- max(1, floor(n_depto_total * prop_muestra))
     }
     
     idx_muestra <- c(
          idx_muestra,
          sample(idx_depto, size = n_depto, replace = FALSE)
     )
}

# Final dataset
dat_muestra <- dat[idx_muestra, , drop = FALSE]

dim(dat_muestra)
table(dat_muestra$estu_depto_reside)

# Bogotá data ------------------------------------------------------------------

dat_bogota <- dat_muestra[
     dat_muestra$estu_depto_reside == "BOGOTÁ",
     ,
     drop = FALSE
]

# Mathematics score
y <- dat_bogota$punt_matematicas
y <- as.numeric(y)
y <- y[!is.na(y)]

# Data summary
summary(y)

# Sample size
(n <- length(y))

# Sufficient statistics
(yb <- mean(y))
(s2 <- var(y))

round(yb, 3)
round(s2, 3)
round(sqrt(s2), 3)

# Gibbs sampler ----------------------------------------------------------------

# Hyperparameters
mu0 <- 50
t20 <- 10^2   # Empirical rule: P(|theta - mu0| < 3 tau0) = 99.7%
s20 <- 10^2
nu0 <- 1

# Data quantities
n  <- length(y)
yb <- mean(y)

# Number of iterations
B <- 10000

# Matrix to store the samples
THETA <- matrix(data = NA_real_, nrow = B, ncol = 2)
colnames(THETA) <- c("theta", "sigma2")

# Initialization from the prior distribution
set.seed(123)

theta <- rnorm(1, mean = mu0, sd = sqrt(t20))
sig2  <- 1 / rgamma(1, shape = nu0 / 2, rate = nu0 * s20 / 2)

# Fixed quantities
nun           <- nu0 + n
paso_progreso <- max(1, floor(B / 10))

# Gibbs sampler iterations
for (b in seq_len(B)) {
     # Update theta
     t2n   <- 1 / (1 / t20 + n / sig2)
     mun   <- t2n * (mu0 / t20 + n * yb / sig2)
     theta <- rnorm(1, mean = mun, sd = sqrt(t2n))
     
     # Update sigma^2
     s2n  <- (nu0 * s20 + sum((y - theta)^2)) / nun
     sig2 <- 1 / rgamma(1, shape = nun / 2, rate = nun * s2n / 2)
     
     # Store samples
     THETA[b, ] <- c(theta, sig2)
     
     # Display progress every 10% of the iterations
     if (b %% paso_progreso == 0) {
          cat(sprintf("%.0f%% completado...\n", 100 * b / B))
     }
}

# Convert the samples to a data frame
THETA <- as.data.frame(THETA)

# Sampler visualization --------------------------------------------------------

# Visualization of the algorithm, for illustration only

# Iterations to visualize
iteraciones <- c(15, 100)
max_iter    <- max(iteraciones)

# Common limits for comparing trajectories
x_lim <- range(THETA$theta[1:max_iter])
y_lim <- range(THETA$sigma2[1:max_iter])

for (m1 in iteraciones) {
     pdf(
          file      = paste0("matematicas_bogota_gibbs_iteraciones_", m1, ".pdf"),
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
          x    = THETA$theta[1:m1],
          y    = THETA$sigma2[1:m1],
          type = "l",
          xlim = x_lim,
          ylim = y_lim,
          lty  = 1,
          col  = "gray",
          xlab = expression(theta),
          ylab = expression(sigma^2)
     )
     
     text(
          x      = THETA$theta[1:m1],
          y      = THETA$sigma2[1:m1],
          labels = seq_len(m1)
     )
     
     dev.off()
}

# Chains -----------------------------------------------------------------------

# Chain plot for theta
pdf(
     file      = "matematicas_bogota_gibbs_cadena_theta.pdf",
     width     = 6,
     height    = 4,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = seq_len(nrow(THETA)),
     y    = THETA$theta,
     type = "p",
     pch  = 16,
     col  = adjustcolor("black", alpha.f = 0.1),
     cex  = 0.5,
     xlab = "Iteration",
     ylab = expression(theta),
     main = ""
)

dev.off()

# Chain plot for sigma^2
pdf(
     file      = "matematicas_bogota_gibbs_cadena_sigma2.pdf",
     width     = 6,
     height    = 4,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = seq_len(nrow(THETA)),
     y    = THETA$sigma2,
     type = "p",
     pch  = 16,
     col  = adjustcolor("black", alpha.f = 0.1),
     cex  = 0.5,
     xlab = "Iteration",
     ylab = expression(sigma^2),
     main = ""
)

dev.off()

# Autocorrelation functions ----------------------------------------------------

# Maximum number of lags
max_lag <- 30

# Autocorrelation function for theta
pdf(
     file      = "matematicas_bogota_gibbs_acf_theta.pdf",
     width     = 6,
     height    = 4,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

acf_theta <- acf(
     x       = THETA$theta,
     lag.max = max_lag,
     main    = "",
     xlab    = "Lag",
     ylab    = "Autocorrelation",
     col     = "black",
     lwd     = 2,
     ci.col  = "gray40"
)

dev.off()

# Autocorrelation function for sigma^2
pdf(
     file      = "matematicas_bogota_gibbs_acf_sigma2.pdf",
     width     = 6,
     height    = 4,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

acf_sigma2 <- acf(
     x       = THETA$sigma2,
     lag.max = max_lag,
     main    = "",
     xlab    = "Lag",
     ylab    = "Autocorrelation",
     col     = "black",
     lwd     = 2,
     ci.col  = "gray40"
)

dev.off()

# Table with the first autocorrelations
tab_acf <- data.frame(
     rezago     = seq_len(10) - 1,
     acf_theta  = as.numeric(acf_theta$acf)[1:10],
     acf_sigma2 = as.numeric(acf_sigma2$acf)[1:10]
)

round(tab_acf, 4)

# Gibbs sampler ----------------------------------------------------------------

mcmc <- function(y, B, mu0, t20, nu0, s20, verbose = TRUE) {
     # Data quantities
     n  <- length(y)
     yb <- mean(y)
     
     # Initialization from the prior distribution
     theta <- rnorm(1, mean = mu0, sd = sqrt(t20))
     sig2  <- 1 / rgamma(1, shape = nu0 / 2, rate = nu0 * s20 / 2)
     
     # Matrix to store the samples
     THETA <- matrix(data = NA_real_, nrow = B, ncol = 2)
     colnames(THETA) <- c("theta", "sigma2")
     
     # Fixed quantities
     nun           <- nu0 + n
     paso_progreso <- max(1, floor(B / 10))
     
     # Gibbs sampler iterations
     for (b in seq_len(B)) {
          # Update theta
          t2n   <- 1 / (1 / t20 + n / sig2)
          mun   <- t2n * (mu0 / t20 + n * yb / sig2)
          theta <- rnorm(1, mean = mun, sd = sqrt(t2n))
          
          # Update sigma^2
          s2n  <- (nu0 * s20 + sum((y - theta)^2)) / nun
          sig2 <- 1 / rgamma(1, shape = nun / 2, rate = nun * s2n / 2)
          
          # Store samples
          THETA[b, ] <- c(theta, sig2)
          
          # Display progress every 10% of the iterations
          if (verbose && b %% paso_progreso == 0) {
               cat(sprintf("%.0f%% completado...\n", 100 * b / B))
          }
     }
     
     as.data.frame(THETA)
}

# Chains -----------------------------------------------------------------------

suppressMessages(suppressWarnings(library(coda)))

set.seed(123)

cadena_1 <- mcmc(y, B, mu0, t20, nu0, s20, verbose = FALSE)
cadena_2 <- mcmc(y, B, mu0, t20, nu0, s20, verbose = FALSE)
cadena_3 <- mcmc(y, B, mu0, t20, nu0, s20, verbose = FALSE)

# Create an mcmc.list object with multiple chains
chains_mcmc <- coda::mcmc.list(
     coda::mcmc(data = as.matrix(cadena_1), start = 1, end = B, thin = 1),
     coda::mcmc(data = as.matrix(cadena_2), start = 1, end = B, thin = 1),
     coda::mcmc(data = as.matrix(cadena_3), start = 1, end = B, thin = 1)
)

# Gelman-Rubin diagnostic
rhat <- coda::gelman.diag(chains_mcmc, autoburnin = FALSE)$psrf

# Organize R-hat results
tab_rhat <- data.frame(
     parametro   = rownames(rhat),
     rhat        = round(rhat[, "Point est."], 4),
     rhat_ic_sup = round(rhat[, "Upper C.I."], 4),
     row.names   = NULL
)

# Numerical diagnostics --------------------------------------------------------

# Required library
suppressMessages(suppressWarnings(library(coda)))

# Convert the samples to mcmc format
THETA_mcmc <- coda::mcmc(THETA)

# Effective sample size
neff <- coda::effectiveSize(THETA_mcmc)

# Monte Carlo standard error
eemc <- apply(X = THETA, MARGIN = 2, FUN = sd) / sqrt(neff)

# Monte Carlo coefficient of variation
cvmc <- 100 * eemc / abs(colMeans(THETA))

# Summary table
tab_diagnosticos <- data.frame(
     parametro = c("theta", "sigma2"),
     neff      = round(as.numeric(neff), 1),
     eemc      = round(as.numeric(eemc), 4),
     cvmc      = round(as.numeric(cvmc), 4)
)

# Merge numerical diagnostics and R-hat
tab_diagnosticos <- merge(
     x     = tab_diagnosticos,
     y     = tab_rhat,
     by    = "parametro",
     all.x = TRUE,
     sort  = FALSE
)

# Organize parameters as rows
rownames(tab_diagnosticos) <- tab_diagnosticos$parametro
tab_diagnosticos$parametro <- NULL

tab_diagnosticos

# Posterior inference ----------------------------------------------------------

# Posterior samples of sigma
sigma <- sqrt(THETA$sigma2)

# Posterior predictive distribution
set.seed(1234)

y_new <- rnorm(
     n    = nrow(THETA),
     mean = THETA$theta,
     sd   = sigma
)

# Summary table
tab <- rbind(
     theta = c(
          mean(THETA$theta),
          sd(THETA$theta) / abs(mean(THETA$theta)),
          quantile(THETA$theta, probs = c(0.025, 0.975))
     ),
     sigma = c(
          mean(sigma),
          sd(sigma) / abs(mean(sigma)),
          quantile(sigma, probs = c(0.025, 0.975))
     ),
     y_pred = c(
          mean(y_new),
          sd(y_new) / abs(mean(y_new)),
          quantile(y_new, probs = c(0.025, 0.975))
     )
)

# Column and row names
colnames(tab) <- c("Estimación", "CV", "L. Inf. 95%", "L. Sup. 95%")
rownames(tab) <- c("Media", "Desviación estándar", "Predicción")

round(tab, 3)

# End --------------------------------------------------------------------------