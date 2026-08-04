# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
dir_work <- "~/Dropbox/UN/bayes_book_en"
setwd(dir_work)

# Libraries
suppressMessages(suppressWarnings(library(mvtnorm)))

# Data ------------------------------------------------------------------------

# Reading comprehension scores before and after instruction
pre_test <- c(
     59, 43, 34, 32, 42, 38, 55, 67, 64, 45, 49,
     72, 34, 70, 34, 50, 41, 52, 60, 34, 28, 35
)

post_test <- c(
     77, 39, 46, 26, 38, 43, 68, 86, 77, 60, 50,
     59, 38, 48, 55, 58, 54, 60, 75, 47, 48, 33
)

# Data matrix
Y <- cbind(
     pre_test  = pre_test,
     post_test = post_test
)

# Dimensions
n <- nrow(Y)
p <- ncol(Y)

# Descriptive inspection
summary(Y)

# Sufficient statistics -------------------------------------------------------

# Sample mean
yb <- colMeans(Y)

# Sample covariance matrix
SS <- cov(Y)

# Results
round(yb, 1)
round(SS, 1)

# Scatter plot -----------------------------------------------------------------

pdf(
     file      = "comprension_lectura_dispersion.pdf",
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
     x    = Y[, "pre_test"],
     y    = Y[, "post_test"],
     pch  = 16,
     col  = adjustcolor(4, alpha.f = 0.7),
     cex  = 1.1,
     xlab = "Before",
     ylab = "After",
     main = ""
)

# Reference lines at the sample means
abline(
     v   = mean(Y[, "pre_test"]),
     col = "gray70",
     lwd = 2,
     lty = 2
)

abline(
     h   = mean(Y[, "post_test"]),
     col = "gray70",
     lwd = 2,
     lty = 2
)

# Point corresponding to the sample mean vector
points(
     x   = mean(Y[, "pre_test"]),
     y   = mean(Y[, "post_test"]),
     pch = 16,
     col = "gray70",
     cex = 1.3
)

dev.off()

# Sample correlation
cor_y <- cor(Y[, "pre_test"], Y[, "post_test"])

round(cor_y, 3)

# Hyperparameters --------------------------------------------------------------

mu0 <- c(50, 50)

L0 <- matrix(
     data  = c(
          278, 139,
          139, 278
     ),
     nrow  = 2,
     ncol  = 2,
     byrow = TRUE
)

nu0 <- 4

S0 <- matrix(
     data  = c(
          278, 139,
          139, 278
     ),
     nrow  = 2,
     ncol  = 2,
     byrow = TRUE
)

# Auxiliary function -----------------------------------------------------------

# Generate a matrix from WI(nu, S^{-1})
riwishart <- function(nu, S) {
     solve(stats::rWishart(n = 1, df = nu, Sigma = solve(S))[, , 1])
}

# Gibbs sampler ----------------------------------------------------------------

mcmc <- function(Y, mu0, L0, nu0, S0, B, seed = 123, verbose = TRUE) {
     # Data quantities
     n  <- nrow(Y)
     p  <- ncol(Y)
     yb <- colMeans(Y)
     SS <- cov(Y)
     
     # Fixed prior quantities
     iL0   <- solve(L0)
     L0mu0 <- iL0 %*% mu0
     nun   <- nu0 + n
     SSn   <- S0 + (n - 1) * SS
     
     # Storage
     THETA <- matrix(data = NA_real_, nrow = B, ncol = p)
     colnames(THETA) <- paste0("theta_", 1:p)
     
     idx_sigma <- expand.grid(row = seq_len(p), col = seq_len(p))
     
     SIGMA <- matrix(data = NA_real_, nrow = B, ncol = p * p)
     colnames(SIGMA) <- paste0("sigma2_", idx_sigma$row, idx_sigma$col)
     
     YS <- matrix(data = NA_real_, nrow = B, ncol = p)
     colnames(YS) <- paste0("y_pred_", colnames(Y))
     
     LL <- numeric(B)
     
     # Initialization from the prior distribution
     Sigma <- riwishart(nu = nu0, S = S0)
     
     # Quantity used to display progress
     paso_progreso <- max(1, floor(B / 10))
     
     # Chain
     set.seed(seed)
     
     for (b in seq_len(B)) {
          # Update theta
          iSigma <- solve(Sigma)
          Ln     <- solve(iL0 + n * iSigma)
          mun    <- Ln %*% (L0mu0 + n * iSigma %*% yb)
          
          theta <- as.numeric(
               mvtnorm::rmvnorm(
                    n     = 1,
                    mean  = as.numeric(mun),
                    sigma = Ln
               )
          )
          
          # Update Sigma
          Sn    <- SSn + n * tcrossprod(yb - theta)
          Sigma <- riwishart(nu = nun, S = Sn)
          
          # Posterior predictive distribution
          YS[b, ] <- as.numeric(
               mvtnorm::rmvnorm(
                    n     = 1,
                    mean  = theta,
                    sigma = Sigma
               )
          )
          
          # Log-likelihood
          LL[b] <- sum(
               mvtnorm::dmvnorm(
                    x     = Y,
                    mean  = theta,
                    sigma = Sigma,
                    log   = TRUE
               )
          )
          
          # Storage
          THETA[b, ] <- theta
          SIGMA[b, ] <- as.vector(Sigma)
          
          if (verbose && b %% paso_progreso == 0) {
               cat(sprintf("%.0f%% completado...\n", 100 * b / B))
          }
     }
     
     list(
          THETA = as.data.frame(THETA),
          SIGMA = as.data.frame(SIGMA),
          LL    = LL,
          YS    = as.data.frame(YS)
     )
}

# Model fitting ----------------------------------------------------------------

fit <- mcmc(
     Y       = Y,
     mu0     = mu0,
     L0      = L0,
     nu0     = nu0,
     S0      = S0,
     B       = 10000,
     seed    = 123,
     verbose = TRUE
)

# Extract results --------------------------------------------------------------

THETA <- fit$THETA
SIGMA <- fit$SIGMA
LL    <- fit$LL
YS    <- fit$YS

# Function for plotting chains -------------------------------------------------

plot_cadena <- function(x, ylab, col, cex, file) {
     pdf(
          file      = file,
          width     = 6,
          height    = 4,
          pointsize = 17
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3.1, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     plot(
          x    = seq_along(x),
          y    = x,
          type = "p",
          pch  = 16,
          col  = adjustcolor(col, alpha.f = 0.1),
          cex  = cex,
          xlab = "Iteration",
          ylab = ylab,
          main = ""
     )
     
     dev.off()
}

# Monte Carlo precision metrics -----------------------------------------------

metricas_mcmc <- function(x) {
     # Convert the chain to mcmc format
     x_mcmc <- coda::mcmc(x)
     
     # Compute metrics
     neff <- round(as.numeric(coda::effectiveSize(x_mcmc)), 1)
     eemc <- round(stats::sd(x) / sqrt(neff), 4)
     cvmc <- round(eemc / abs(mean(x)), 4)
     
     c(
          neff = neff,
          eemc = eemc,
          cvmc = cvmc
     )
}

# R-hat diagnostic with multiple chains ---------------------------------------

calcular_rhat <- function(Y, mu0, L0, nu0, S0, B, M = 3, seed = 123) {
     # Run multiple chains
     ajustes <- vector(mode = "list", length = M)
     
     for (m in seq_len(M)) {
          ajustes[[m]] <- mcmc(
               Y       = Y,
               mu0     = mu0,
               L0      = L0,
               nu0     = nu0,
               S0      = S0,
               B       = B,
               seed    = seed + m,
               verbose = FALSE
          )
     }
     
     # Convert THETA to mcmc.list
     cadenas_theta <- coda::mcmc.list(
          lapply(
               X   = ajustes,
               FUN = function(x) coda::mcmc(as.matrix(x$THETA))
          )
     )
     
     # Convert SIGMA to mcmc.list
     cadenas_sigma <- coda::mcmc.list(
          lapply(
               X   = ajustes,
               FUN = function(x) coda::mcmc(as.matrix(x$SIGMA))
          )
     )
     
     # Convert LL to mcmc.list
     cadenas_ll <- coda::mcmc.list(
          lapply(
               X   = ajustes,
               FUN = function(x) {
                    coda::mcmc(matrix(x$LL, ncol = 1, dimnames = list(NULL, "LL")))
               }
          )
     )
     
     # Compute R-hat
     rhat_theta <- coda::gelman.diag(
          x          = cadenas_theta,
          autoburnin = FALSE
     )$psrf
     
     rhat_sigma <- coda::gelman.diag(
          x          = cadenas_sigma,
          autoburnin = FALSE
     )$psrf
     
     rhat_ll <- coda::gelman.diag(
          x          = cadenas_ll,
          autoburnin = FALSE
     )$psrf
     
     # Organize results
     tab_theta <- data.frame(
          parametro   = rownames(rhat_theta),
          rhat        = rhat_theta[, "Point est."],
          rhat_ic_sup = rhat_theta[, "Upper C.I."],
          tipo        = "theta",
          row.names   = NULL
     )
     
     tab_sigma <- data.frame(
          parametro   = rownames(rhat_sigma),
          rhat        = rhat_sigma[, "Point est."],
          rhat_ic_sup = rhat_sigma[, "Upper C.I."],
          tipo        = "sigma",
          row.names   = NULL
     )
     
     tab_ll <- data.frame(
          parametro   = rownames(rhat_ll),
          rhat        = rhat_ll[, "Point est."],
          rhat_ic_sup = rhat_ll[, "Upper C.I."],
          tipo        = "logverosimilitud",
          row.names   = NULL
     )
     
     tab_rhat <- rbind(tab_theta, tab_sigma, tab_ll)
     
     rownames(tab_rhat) <- tab_rhat$parametro
     tab_rhat$parametro <- NULL
     
     list(
          rhat    = tab_rhat,
          ajustes = ajustes
     )
}

# Posterior summary ------------------------------------------------------------

resumen_posterior <- function(x, nivel_confianza = 0.95) {
     # Probabilities for the credible interval
     alpha <- 1 - nivel_confianza
     probs <- c(alpha / 2, 1 - alpha / 2)
     
     # Posterior summary
     media <- mean(x)
     cv    <- sd(x) / abs(media)
     ic    <- quantile(x, probs = probs, names = FALSE)
     
     round(
          c(
               media = media,
               cv    = cv,
               li    = ic[1],
               ls    = ic[2]
          ),
          digits = 3
     )
}

# Convergence ------------------------------------------------------------------

# Chains for the log-likelihood
plot_cadena(
     x    = LL,
     ylab = "Log-likelihood",
     col  = 1,
     cex  = 0.5,
     file = "comprension_lectura_cadena_logverosimilitud.pdf"
)

# Chains for each component of theta
ylab_theta <- list(
     expression(theta[1]),
     expression(theta[2])
)

for (j in seq_along(THETA)) {
     plot_cadena(
          x    = THETA[[j]],
          ylab = ylab_theta[[j]],
          col  = 1,
          cex  = 0.5,
          file = paste0("comprension_lectura_cadena_", colnames(THETA)[j], ".pdf")
     )
}

# Chains for each component of Sigma
ylab_sigma <- list(
     expression(sigma[1]^2),
     expression(sigma[21]),
     expression(sigma[12]),
     expression(sigma[2]^2)
)

for (j in seq_along(SIGMA)) {
     plot_cadena(
          x    = SIGMA[[j]],
          ylab = ylab_sigma[[j]],
          col  = 1,
          cex  = 0.5,
          file = paste0("comprension_lectura_cadena_", colnames(SIGMA)[j], ".pdf")
     )
}

# Convergence metrics for the log-likelihood
metricas_mcmc(LL)

# Convergence metrics for THETA
tab_theta <- t(
     sapply(
          X   = THETA,
          FUN = metricas_mcmc
     )
)

round(tab_theta, 5)

# Convergence metrics for SIGMA
tab_sigma <- t(
     sapply(
          X   = SIGMA,
          FUN = metricas_mcmc
     )
)

round(tab_sigma, 5)

# R-hat
diagnostico_rhat <- calcular_rhat(
     Y    = Y,
     mu0  = mu0,
     L0   = L0,
     nu0  = nu0,
     S0   = S0,
     B    = 10000,
     M    = 3,
     seed = 123
)

round(diagnostico_rhat$rhat[, c("rhat", "rhat_ic_sup")], 4)

# Inference --------------------------------------------------------------------

# Posterior difference between means
delta_theta <- THETA$theta_2 - THETA$theta_1

# Posterior predictive difference
delta_pred <- YS$y_pred_post_test - YS$y_pred_pre_test

# Posterior correlation induced by Sigma
rho_sigma <- SIGMA$sigma2_12 / sqrt(SIGMA$sigma2_11 * SIGMA$sigma2_22)

# Posterior summaries
tab_resumen <- rbind(
     delta_theta = resumen_posterior(
          x               = delta_theta,
          nivel_confianza = 0.95
     ),
     delta_pred = resumen_posterior(
          x               = delta_pred,
          nivel_confianza = 0.95
     ),
     rho_sigma = resumen_posterior(
          x               = rho_sigma,
          nivel_confianza = 0.95
     )
)

tab_resumen

# Posterior probabilities
prob_delta_theta_pos <- mean(delta_theta > 0)
prob_delta_pred_pos  <- mean(delta_pred > 0)

round(prob_delta_theta_pos, 4)
round(prob_delta_pred_pos, 4)

# Posterior plot of theta_2 vs theta_1 -----------------------------------------

lim_theta <- range(THETA[, 1], THETA[, 2])

pdf(
     file      = "comprension_lectura_posterior_delta_theta.pdf",
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
     x    = THETA[, 1],
     y    = THETA[, 2],
     pch  = 16,
     col  = adjustcolor(4, alpha.f = 0.1),
     cex  = 0.5,
     xlim = lim_theta,
     ylim = lim_theta,
     xlab = expression(theta[1]),
     ylab = expression(theta[2]),
     main = "",
     asp  = 1
)

abline(
     v   = mean(THETA[, 1]),
     col = "gray75",
     lwd = 2,
     lty = 2
)

abline(
     h   = mean(THETA[, 2]),
     col = "gray75",
     lwd = 2,
     lty = 2
)

abline(
     a   = 0,
     b   = 1,
     col = "gray30",
     lwd = 2
)

points(
     x   = mean(THETA[, 1]),
     y   = mean(THETA[, 2]),
     pch = 3,
     col = 2,
     lwd = 2,
     cex = 1.3
)

dev.off()

# Posterior predictive plot of y_2^* vs y_1^* ---------------------------------

lim_pred <- range(YS[, 1], YS[, 2])

pdf(
     file      = "comprension_lectura_posterior_delta_pred.pdf",
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
     x    = YS[, 1],
     y    = YS[, 2],
     pch  = 16,
     col  = adjustcolor(4, alpha.f = 0.1),
     cex  = 0.5,
     xlim = lim_pred,
     ylim = lim_pred,
     xlab = expression(tilde(y)[1]),
     ylab = expression(tilde(y)[2]),
     main = "",
     asp  = 1
)

abline(
     v   = mean(YS[, 1]),
     col = "gray75",
     lwd = 2,
     lty = 2
)

abline(
     h   = mean(YS[, 2]),
     col = "gray75",
     lwd = 2,
     lty = 2
)

abline(
     a   = 0,
     b   = 1,
     col = "gray30",
     lwd = 2
)

points(
     x   = mean(YS[, 1]),
     y   = mean(YS[, 2]),
     pch = 3,
     col = 2,
     lwd = 2,
     cex = 1.3
)

dev.off()

# Posterior distribution of the correlation -----------------------------------

RHO <- SIGMA[, 3] / sqrt(SIGMA[, 1] * SIGMA[, 4])

pdf(
     file      = "comprension_lectura_posterior_correlacion.pdf",
     width     = 5,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

hist(
     x      = RHO,
     freq   = FALSE,
     breaks = 30,
     col    = "gray90",
     border = "white",
     xlim   = c(0, 1),
     xlab   = expression(rho),
     ylab   = "Density",
     main   = ""
)

abline(
     v   = mean(RHO),
     col = 2,
     lwd = 2,
     lty = 2
)

abline(
     v   = quantile(RHO, probs = c(0.025, 0.975)),
     col = 4,
     lwd = 2,
     lty = 2
)

legend(
     x      = "topleft",
     legend = c("Mean", "95% CI"),
     col    = c(2, 4),
     fill   = c(2, 4),
     border = c(2, 4),
     bty    = "n"
)

dev.off()

# End --------------------------------------------------------------------------