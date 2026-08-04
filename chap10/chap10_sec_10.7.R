# Settings ---------------------------------------------------------------------

rm(list = ls())

setwd("~/Dropbox/UN/bayes_book_en")

# Data -------------------------------------------------------------------------

# Dataset
spfage <- structure(
     c(
          3, 1, 1, 2, 0, 0, 6, 3, 4, 2, 1, 6, 2, 3, 3, 4, 7, 2, 2, 1,
          1, 3, 5, 5, 0, 2, 1, 2, 6, 6, 2, 2, 0, 2, 4, 1, 2, 5, 1, 2,
          1, 0, 0, 2, 4, 2, 2, 2, 2, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
          2, 2, 2, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4,
          4, 4, 5, 5, 5, 5, 3, 3, 3, 3, 3, 3, 3, 6, 1, 1, 9, 9, 1, 1,
          1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 25, 25, 16, 16,
          16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 25, 16, 16, 16, 16,
          25, 25, 25, 25, 9, 9, 9, 9, 9, 9, 9, 36, 1, 1
     ),
     .Dim = c(52L, 4L),
     .Dimnames = list(
          NULL,
          c("fledged", "intercept", "age", "age2")
     )
)

spfage <- as.data.frame(spfage)

y <- spfage$fledged                     # response variable
X <- cbind(1, spfage$age, spfage$age2)  # design matrix
n <- dim(X)[1]                          # sample size
p <- dim(X)[2]

# Descriptive plots ------------------------------------------------------------

tabla_crias <- data.frame(
     crias      = as.integer(names(table(spfage$fledged))),
     frecuencia = as.integer(table(spfage$fledged)),
     proporcion = as.numeric(prop.table(table(spfage$fledged))),
     row.names  = NULL
)

tabla_edad <- data.frame(
     edad       = as.integer(names(table(spfage$age))),
     frecuencia = as.integer(table(spfage$age)),
     proporcion = as.numeric(prop.table(table(spfage$age))),
     row.names  = NULL
)

resumen_por_edad <- do.call(
     data.frame,
     aggregate(
          fledged ~ age,
          data = spfage,
          FUN  = function(x) {
               c(
                    n       = length(x),
                    minimo  = min(x),
                    q1      = quantile(x, 0.25),
                    mediana = median(x),
                    media   = mean(x),
                    var     = var(x),
                    sd      = sd(x),
                    q3      = quantile(x, 0.75),
                    maximo  = max(x)
               )
          }
     )
)

names(resumen_por_edad) <- c(
     "edad",
     "n",
     "minimo",
     "q1",
     "mediana",
     "media",
     "var",
     "sd",
     "q3",
     "maximo"
)

# Distribution of the number of fledglings

pdf(
     file      = "gorriones_distribucion_crias.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

barplot(
     height    = tabla_crias$frecuencia,
     names.arg = tabla_crias$crias,
     border    = "white",
     col       = "gray",
     xlab      = "Number of fledglings",
     ylab      = "Frequency",
     main      = ""
)

box()

dev.off()

# Age distribution

pdf(
     file      = "gorriones_distribucion_edad.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

barplot(
     height    = tabla_edad$frecuencia,
     names.arg = tabla_edad$edad,
     border    = "white",
     col       = "gray",
     xlab      = "Age (years)",
     ylab      = "Frequency",
     main      = ""
)

box()

dev.off()

# Scatter plot

pdf(
     file      = "gorriones_dispersion_crias_edad.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

plot(
     x    = jitter(spfage$age, amount = 0.06),
     y    = spfage$fledged,
     pch  = 16,
     col  = adjustcolor(1, 0.6),
     xlab = "Age (years)",
     ylab = "Number of fledglings",
     main = "",
     xaxt = "n"
)

axis(
     side = 1,
     at   = sort(unique(spfage$age))
)

lines(
     x    = resumen_por_edad$edad,
     y    = resumen_por_edad$media,
     type = "b",
     pch  = 16,
     col  = 2,
     lwd  = 2
)

box()

dev.off()

# Boxplots by age

pdf(
     file      = "gorriones_cajas_crias_edad.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

boxplot(
     fledged ~ age,
     data   = spfage,
     xlab   = "Age (years)",
     ylab   = "Number of fledglings",
     col    = "gray",
     border = 1,
     main   = ""
)

points(
     x   = seq_along(resumen_por_edad$edad),
     y   = resumen_por_edad$media,
     pch = 16,
     col = 1
)

box()

dev.off()

# Prior distribution -----------------------------------------------------------

mu0 <- rep(0, p)

Sigma0 <- 10 * diag(p)

Sigma0_inv <- solve(Sigma0)

# Metropolis -------------------------------------------------------------------

# Unnormalized log-posterior function

log_posterior <- function(beta, y, X, mu0, Sigma0_inv) {
     eta <- drop(X %*% beta)
     
     log_likelihood <- sum(
          dpois(
               x      = y,
               lambda = exp(eta),
               log    = TRUE
          )
     )
     
     beta_centered <- beta - mu0
     
     log_prior <- -0.5 * drop(
          t(beta_centered) %*% Sigma0_inv %*% beta_centered
     )
     
     log_likelihood + log_prior
}

# Algorithm configuration

S <- 11000

Delta <- 0.7 * solve(crossprod(X))

beta_current <- rep(0, p)

log_post_current <- log_posterior(
     beta       = beta_current,
     y          = y,
     X          = X,
     mu0        = mu0,
     Sigma0_inv = Sigma0_inv
)

BET_metropolis <- matrix(
     data = NA,
     nrow = S,
     ncol = p
)

colnames(BET_metropolis) <- paste0("beta", seq_len(p))

accept <- 0

# Chain

set.seed(123)

for (s in seq_len(S)) {
     beta_proposed <- drop(
          mvtnorm::rmvnorm(
               n     = 1,
               mean  = beta_current,
               sigma = Delta
          )
     )
     
     log_post_proposed <- log_posterior(
          beta       = beta_proposed,
          y          = y,
          X          = X,
          mu0        = mu0,
          Sigma0_inv = Sigma0_inv
     )
     
     log_r <- log_post_proposed - log_post_current
     
     if (log(runif(1)) <= log_r) {
          beta_current     <- beta_proposed
          log_post_current <- log_post_proposed
          accept           <- accept + 1
     }
     
     BET_metropolis[s, ] <- beta_current
     
     if (s %% floor(S / 10) == 0) {
          cat(100 * round(s / S, 1), "%\n", sep = "")
     }
}

# Acceptance rate

acceptance_rate <- accept / S

round(acceptance_rate, 3)

# HMC --------------------------------------------------------------------------

# Gradient of the unnormalized log-posterior

grad_log_posterior <- function(beta, y, X, mu0, Sigma0_inv) {
     eta   <- drop(X %*% beta)
     theta <- exp(eta)
     
     grad_log_likelihood <- drop(
          crossprod(X, y - theta)
     )
     
     grad_log_prior <- -drop(
          Sigma0_inv %*% (beta - mu0)
     )
     
     grad_log_likelihood + grad_log_prior
}

# Algorithm configuration

S <- 11000

L <- 100

eps <- 1 / L

beta_current <- rep(0, p)

log_post_current <- log_posterior(
     beta       = beta_current,
     y          = y,
     X          = X,
     mu0        = mu0,
     Sigma0_inv = Sigma0_inv
)

BET_hmc <- matrix(
     data = NA,
     nrow = S,
     ncol = p
)

colnames(BET_hmc) <- paste0("beta", seq_len(p))

accept <- 0

# HMC algorithm with identity mass matrix

set.seed(123)

for (s in seq_len(S)) {
     beta_proposed <- beta_current
     
     phi_current <- rnorm(
          n    = p,
          mean = 0,
          sd   = 1
     )
     
     phi_proposed <- phi_current
     
     # First half-step for the momentum
     
     phi_proposed <- phi_proposed +
          (eps / 2) * grad_log_posterior(
               beta       = beta_proposed,
               y          = y,
               X          = X,
               mu0        = mu0,
               Sigma0_inv = Sigma0_inv
          )
     
     # Leapfrog steps
     
     for (l in seq_len(L)) {
          beta_proposed <- beta_proposed + eps * phi_proposed
          
          if (l < L) {
               phi_proposed <- phi_proposed +
                    eps * grad_log_posterior(
                         beta       = beta_proposed,
                         y          = y,
                         X          = X,
                         mu0        = mu0,
                         Sigma0_inv = Sigma0_inv
                    )
          }
     }
     
     # Final half-step for the momentum
     
     phi_proposed <- phi_proposed +
          (eps / 2) * grad_log_posterior(
               beta       = beta_proposed,
               y          = y,
               X          = X,
               mu0        = mu0,
               Sigma0_inv = Sigma0_inv
          )
     
     # Reverse the momentum to ensure reversibility
     
     phi_proposed <- -phi_proposed
     
     log_post_proposed <- log_posterior(
          beta       = beta_proposed,
          y          = y,
          X          = X,
          mu0        = mu0,
          Sigma0_inv = Sigma0_inv
     )
     
     log_joint_current <- log_post_current +
          sum(dnorm(phi_current, log = TRUE))
     
     log_joint_proposed <- log_post_proposed +
          sum(dnorm(phi_proposed, log = TRUE))
     
     log_r <- log_joint_proposed - log_joint_current
     
     if (log(runif(1)) <= log_r) {
          beta_current     <- beta_proposed
          log_post_current <- log_post_proposed
          accept           <- accept + 1
     }
     
     BET_hmc[s, ] <- beta_current
     
     if (s %% floor(S / 10) == 0) {
          cat(100 * round(s / S, 1), "%\n", sep = "")
     }
}

# Acceptance rate

acceptance_rate_hmc <- accept / S

round(acceptance_rate_hmc, 3)

# Posterior samples ------------------------------------------------------------

burn_in <- 1000

BET_metropolis_post <- BET_metropolis[-seq_len(burn_in), , drop = FALSE]
BET_hmc_post        <- BET_hmc[-seq_len(burn_in), , drop = FALSE]

S_post <- nrow(BET_metropolis_post)

# Chain plots ------------------------------------------------------------------

for (j in seq_len(p)) {
     rango_j <- range(
          BET_metropolis_post[, j],
          BET_hmc_post[, j]
     )
     
     pdf(
          file      = paste0("gorriones_metropolis_cadena_beta_", j, ".pdf"),
          width     = 5.5,
          height    = 3.5,
          pointsize = 15
     )
     
     par(
          mar = c(3, 3, 1.4, 1.4),
          mgp = c(1.75, 0.75, 0)
     )
     
     plot(
          x        = seq_len(S_post),
          y        = BET_metropolis_post[, j],
          type     = "p",
          pch      = 16,
          cex      = 0.3,
          col      = adjustcolor(1, 0.7),
          cex.axis = 0.8,
          ylim     = rango_j,
          xlab     = "Iteration",
          ylab     = bquote(beta[.(j)]),
          main     = ""
     )
     
     box()
     
     dev.off()
     
     pdf(
          file      = paste0("gorriones_hmc_cadena_beta_", j, ".pdf"),
          width     = 5.5,
          height    = 3.5,
          pointsize = 15
     )
     
     par(
          mar = c(3, 3, 1.4, 1.4),
          mgp = c(1.75, 0.75, 0)
     )
     
     plot(
          x        = seq_len(S_post),
          y        = BET_hmc_post[, j],
          type     = "p",
          pch      = 16,
          cex      = 0.3,
          col      = adjustcolor(4, 0.7),
          cex.axis = 0.8,
          ylim     = rango_j,
          xlab     = "Iteration",
          ylab     = bquote(beta[.(j)]),
          main     = ""
     )
     
     box()
     
     dev.off()
}

# Effective sample sizes and Monte Carlo standard errors -----------------------

ess_metropolis <- coda::effectiveSize(BET_metropolis_post)
ess_hmc        <- coda::effectiveSize(BET_hmc_post)

mcse_metropolis <- apply(BET_metropolis_post, 2, sd) / sqrt(ess_metropolis)
mcse_hmc        <- apply(BET_hmc_post, 2, sd) / sqrt(ess_hmc)

tabla_diagnosticos <- data.frame(
     Parametro        = paste0("beta_", seq_len(p)),
     ESS_Metropolis   = round(ess_metropolis, 1),
     ESS_HMC          = round(ess_hmc, 1),
     MCSE_Metropolis  = round(mcse_metropolis, 3),
     MCSE_HMC         = round(mcse_hmc, 3),
     row.names        = NULL
)

tabla_diagnosticos

# Inference --------------------------------------------------------------------

est <- colMeans(x = BET_hmc_post)

sd_post <- apply(
     X      = BET_hmc_post,
     MARGIN = 2,
     FUN    = sd
)

ci <- apply(
     X      = BET_hmc_post,
     MARGIN = 2,
     FUN    = function(x) quantile(x, c(0.025, 0.975))
)

tab <- rbind(
     media = est,
     sd    = sd_post,
     ci
)

round(t(tab), 3)

# Posterior distribution of beta_1 ---------------------------------------------

pdf(
     file      = paste0("gorriones_post_beta_1.pdf"),
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

hist(
     BET_hmc_post[, 1],
     freq   = FALSE,
     border = "white",
     col    = "gray",
     main   = "",
     nclass = 20,
     xlim   = c(-1, 2),
     ylim   = c(0, 7),
     xlab   = expression(beta[1]),
     ylab   = "Density"
)

abline(v = est[1], lty = 1, col = 2)
abline(v = ci[, 1], lty = 2, col = 4)

dev.off()

# Posterior distribution of beta_2 ---------------------------------------------

pdf(
     file      = paste0("gorriones_post_beta_2.pdf"),
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

hist(
     BET_hmc_post[, 2],
     freq   = FALSE,
     border = "white",
     col    = "gray",
     main   = "",
     nclass = 20,
     xlim   = c(-1, 2),
     ylim   = c(0, 7),
     xlab   = expression(beta[2]),
     ylab   = "Density"
)

abline(v = est[2], lty = 1, col = 2)
abline(v = ci[, 2], lty = 2, col = 4)

dev.off()

# Posterior distribution of beta_3 ---------------------------------------------

pdf(
     file      = paste0("gorriones_post_beta_3.pdf"),
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

hist(
     BET_hmc_post[, 3],
     freq   = FALSE,
     border = "gray",
     col    = "gray",
     main   = "",
     nclass = 20,
     xlim   = c(-1, 2),
     ylim   = c(0, 7),
     xlab   = expression(beta[3]),
     ylab   = "Density"
)

abline(v = est[3], lty = 1, col = 2)
abline(v = ci[, 3], lty = 2, col = 4)

dev.off()

# Expected number of fledglings by age -----------------------------------------

pdf(
     file      = paste0("gorriones_post_crias_vs_edad.pdf"),
     width     = 5,
     height    = 5,
     pointsize = 15
)

XX <- cbind(
     rep(1, 6),
     1:6,
     (1:6)^2
)

eXB <- exp(t(XX %*% t(BET_hmc_post)))

qE <- apply(
     X      = eXB,
     MARGIN = 2,
     FUN    = quantile,
     probs  = c(0.025, 0.5, 0.975)
)

par(
     mar = c(3, 3, 1.4, 1.4),
     mgp = c(1.75, 0.75, 0)
)

plot(
     x    = 1:6,
     y    = qE[2, ],
     type = "p",
     pch  = 16,
     ylim = range(qE),
     xlab = "Age (years)",
     ylab = "Number of fledglings"
)

for (j in 1:6) {
     segments(
          x0 = j,
          y0 = qE[1, j],
          x1 = j,
          y1 = qE[3, j]
     )
}

dev.off()

# End --------------------------------------------------------------------------