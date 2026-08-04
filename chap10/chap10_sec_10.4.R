# Settings ---------------------------------------------------------------------

rm(list = ls())

setwd("~/Dropbox/UN/bayes_book_en")

# Libraries
library(readxl)
library(dplyr)
library(stringr)
library(stringi)
library(mvtnorm)
library(coda)

# Data -------------------------------------------------------------------------

# Import files
ideal_points <- read_excel("parapolitica_ideal_points.xlsx")
respuesta    <- read_excel("parapolitica_respuesta.xlsx")

# Standardize the name of the first column
names(ideal_points)[1] <- "legislador"
names(respuesta)[1]    <- "legislador"

# Function for constructing a matching key
normalizar_nombre <- function(x) {
     x |>
          str_replace("\\s*:.*$", "") |>                 # removes party after :
          str_squish() |>                                # removes extra spaces
          str_to_upper() |>                              # converts to uppercase
          stringi::stri_trans_general("Latin-ASCII") |>  # removes accents and ñ
          str_replace_all("[^A-Z ]", "") |>              # removes unusual symbols
          str_squish()
}

# Create keys
ideal_points <- ideal_points |>
     mutate(id_legislador = normalizar_nombre(legislador))

respuesta <- respuesta |>
     mutate(id_legislador = normalizar_nombre(legislador))

# Review names without a match
sin_match <- ideal_points |>
     anti_join(respuesta, by = "id_legislador") |>
     select(legislador, id_legislador)

sin_match

# Merge datasets
parapolitica_modelo <- ideal_points |>
     left_join(
          respuesta |>
               select(-legislador),
          by = "id_legislador"
     ) |>
     select(
          legislador,
          Mean,
          parapolitica0,
          parapolitica1
     )

# Adjust variable names
parapolitica_modelo <- parapolitica_modelo |>
     rename(
          ideal   = Mean,
          parapol = parapolitica0
     ) |>
     select(
          legislador,
          ideal,
          parapol
     )

# Ensure the correct format for the response variable
parapolitica_modelo <- parapolitica_modelo |>
     mutate(
          parapol = as.integer(parapol)
     )

# Review the result
glimpse(parapolitica_modelo)

# Fit the frequentist logit model ----------------------------------------------

mod_logit_freq <- glm(
     parapol ~ ideal,
     data   = parapolitica_modelo,
     family = binomial(link = "logit")
)

# Model summary
summary(mod_logit_freq)

# Estimated coefficients
coef(mod_logit_freq)

# Odds ratios
exp(coef(mod_logit_freq))

# 95% confidence intervals for the coefficients
confint(mod_logit_freq)

# 95% confidence intervals for the odds ratios
exp(confint(mod_logit_freq))

# Data for model fitting -------------------------------------------------------

# Response variable
y <- parapolitica_modelo$parapol

# Design matrix: intercept and ideal point
X <- cbind(
     1,
     ideal = parapolitica_modelo$ideal
)

# Range of the ideal points
round(range(parapolitica_modelo$ideal), 3)

# Sample size and number of predictors
n <- nrow(X)
p <- ncol(X)

# Exploratory analysis ---------------------------------------------------------

# Distribution of the response variable
tabla_y <- table(y)

prop_y <- prop.table(tabla_y)

resumen_y <- data.frame(
     parapolitica = names(tabla_y),
     frecuencia   = as.numeric(tabla_y),
     proporcion   = round(as.numeric(prop_y), 3)
)

resumen_y

# Overall summary of the ideal point
resumen_ideal <- c(
     media   = mean(parapolitica_modelo$ideal),
     mediana = median(parapolitica_modelo$ideal),
     sd      = sd(parapolitica_modelo$ideal),
     minimo  = min(parapolitica_modelo$ideal),
     maximo  = max(parapolitica_modelo$ideal)
)

round(resumen_ideal, 3)

# Summary of the ideal point by parapolitics status
resumen_por_grupo <- aggregate(
     ideal ~ parapol,
     data = parapolitica_modelo,
     FUN  = function(x) {
          c(
               n       = length(x),
               media   = mean(x),
               mediana = median(x),
               sd      = sd(x),
               minimo  = min(x),
               maximo  = max(x)
          )
     }
)

resumen_por_grupo <- do.call(
     data.frame,
     resumen_por_grupo
)

names(resumen_por_grupo) <- c(
     "parapolitica",
     "n",
     "media",
     "mediana",
     "sd",
     "minimo",
     "maximo"
)

round(resumen_por_grupo, 3)

# Histogram of ideal points
pdf(
     file      = "parapolitica_eda_histograma_puntos_ideales.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

hist(
     x      = X[, "ideal"],
     freq   = FALSE,
     col    = "gray85",
     border = "white",
     xlab   = "Ideal point",
     ylab   = "Density",
     main   = ""
)

lines(
     density(X[, "ideal"]),
     lwd = 2
)

rug(
     X[, "ideal"],
     col = adjustcolor("black", 0.4)
)

abline(
     v   = mean(X[, "ideal"]),
     col = 2,
     lwd = 2,
     lty = 2
)

box()

dev.off()

# Boxplot of ideal points by group
pdf(
     file      = "parapolitica_eda_boxplot_puntos_ideales.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

grupo_parapol <- factor(
     parapolitica_modelo$parapol,
     levels = c(0, 1),
     labels = c("No involucrado", "Involucrado")
)

ylim_box <- range(parapolitica_modelo$ideal)
ylim_box <- ylim_box + c(-0.08, 0.08) * diff(ylim_box)

boxplot(
     parapolitica_modelo$ideal ~ grupo_parapol,
     outline  = FALSE,
     col      = "gray90",
     border   = "gray30",
     boxwex   = 0.45,
     lwd      = 1.4,
     medcol   = 1,
     medlwd   = 2,
     whisklty = 1,
     staplelty = 1,
     cex.axis = 0.85,
     xlab     = "",
     ylab     = "Ideal point",
     ylim     = ylim_box,
     main     = ""
)

stripchart(
     parapolitica_modelo$ideal ~ grupo_parapol,
     vertical = TRUE,
     method   = "jitter",
     jitter   = 0.12,
     pch      = 16,
     cex      = 0.55,
     col      = adjustcolor("black", 0.45),
     add      = TRUE
)

points(
     x   = 1:2,
     y   = tapply(parapolitica_modelo$ideal, grupo_parapol, mean),
     pch = 18,
     cex = 1.25,
     col = 1
)

abline(
     h   = 0,
     lty = 3,
     col = "gray60"
)

box()

dev.off()

# MCMC for the logit model -----------------------------------------------------

log1pexp <- function(x) {
     ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

log_likelihood_logit <- function(beta, y, X) {
     eta <- as.numeric(X %*% beta)
     sum(y * eta - log1pexp(eta))
}

log_posterior_logit <- function(beta, y, X, beta0, Sigma0) {
     log_lik <- log_likelihood_logit(beta, y, X)
     
     log_prior <- mvtnorm::dmvnorm(
          x     = beta,
          mean  = beta0,
          sigma = Sigma0,
          log   = TRUE
     )
     
     log_lik + log_prior
}

mcmc_logit <- function(
          y,
          X,
          beta0,
          Sigma0,
          B,
          burn_in,
          thin,
          seed = 123,
          progress = TRUE
) {
     # Basic checks
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     n <- nrow(X)
     p <- ncol(X)
     
     beta_ini <- rep(0, p)
     
     # Proposal calibration
     X_cov <- X[, -1, drop = FALSE]
     
     mod_logit_freq <- glm(
          y ~ X_cov,
          family = binomial(link = "logit")
     )
     
     Delta <- (2.38^2 / p) * vcov(mod_logit_freq)
     
     # Storage
     B_post <- burn_in + B * thin
     
     BETA_post <- matrix(NA, nrow = B, ncol = p)
     LL_post   <- numeric(B)
     
     colnames(BETA_post) <- paste0("beta", seq_len(p) - 1)
     
     ncat <- max(1, floor(B_post / 10))
     
     # Chain
     beta <- beta_ini
     acr  <- 0
     s    <- 0
     
     set.seed(seed)
     
     for (b in seq_len(B_post)) {
          # 1. Proposal
          beta_p <- c(
               mvtnorm::rmvnorm(
                    n     = 1,
                    mean  = beta,
                    sigma = Delta
               )
          )
          
          # 2. Acceptance ratio on the logarithmic scale
          log_post_star <- log_posterior_logit(beta_p, y, X, beta0, Sigma0)
          log_post_curr <- log_posterior_logit(beta, y, X, beta0, Sigma0)
          
          log_r <- log_post_star - log_post_curr
          
          # 3. Accept or reject
          if (log(runif(1)) <= min(0, log_r)) {
               beta <- beta_p
               acr  <- acr + 1
          }
          
          # 4. Store only after warm-up and every thin iterations
          if (b > burn_in && (b - burn_in) %% thin == 0) {
               s <- s + 1
               
               BETA_post[s, ] <- beta
               LL_post[s]     <- log_likelihood_logit(beta, y, X)
          }
          
          # 5. Progress
          if (progress && b %% ncat == 0) {
               cat(round(100 * b / B_post, 1), "% completado ...\n", sep = "")
          }
     }
     
     # Empirical acceptance rate
     tasa_aceptacion <- acr / B_post
     
     # Return results
     list(
          BETA            = BETA_post,
          LL              = LL_post,
          tasa_aceptacion = tasa_aceptacion,
          Delta           = Delta
     )
}

# MCMC for the augmented probit model -----------------------------------------

log_likelihood_probit <- function(beta, y, X) {
     eta <- as.numeric(X %*% beta)
     
     sum(
          ifelse(
               y == 1,
               pnorm(eta, log.p = TRUE),
               pnorm(eta, lower.tail = FALSE, log.p = TRUE)
          )
     )
}

rnorm_trunc_probit <- function(mu, y) {
     n <- length(mu)
     z <- numeric(n)
     
     # y = 1: Normal truncated to (0, Inf)
     id1 <- which(y == 1)
     
     if (length(id1) > 0) {
          z[id1] <- truncnorm::rtruncnorm(
               n    = length(id1),
               a    = 0,
               b    = Inf,
               mean = mu[id1],
               sd   = 1
          )
     }
     
     # y = 0: Normal truncated to (-Inf, 0]
     id0 <- which(y == 0)
     
     if (length(id0) > 0) {
          z[id0] <- truncnorm::rtruncnorm(
               n    = length(id0),
               a    = -Inf,
               b    = 0,
               mean = mu[id0],
               sd   = 1
          )
     }
     
     z
}

mcmc_probit <- function(
          y,
          X,
          beta0,
          Sigma0,
          B,
          burn_in,
          thin,
          seed = 123,
          progress = TRUE
) {
     # Basic checks
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     n <- nrow(X)
     p <- ncol(X)
     
     beta_ini <- rep(0, p)
     
     # Prior precision matrix
     P0 <- solve(Sigma0)
     
     # Conditional covariance matrix of beta
     V_beta <- solve(P0 + t(X) %*% X)
     
     # Storage
     B_post <- burn_in + B * thin
     
     BETA_post <- matrix(NA, nrow = B, ncol = p)
     LL_post   <- numeric(B)
     
     colnames(BETA_post) <- paste0("beta", seq_len(p) - 1)
     
     ncat <- max(1, floor(B_post / 10))
     
     # Chain
     beta <- beta_ini
     z    <- ifelse(y == 1, 1, -1)
     s    <- 0
     
     set.seed(seed)
     
     for (b in seq_len(B_post)) {
          # 1. Simulate latent variables z
          eta <- as.numeric(X %*% beta)
          z   <- rnorm_trunc_probit(mu = eta, y = y)
          
          # 2. Simulate beta from its multivariate Normal full conditional distribution
          m_beta <- V_beta %*% (P0 %*% beta0 + t(X) %*% z)
          
          beta <- c(
               mvtnorm::rmvnorm(
                    n     = 1,
                    mean  = c(m_beta),
                    sigma = V_beta
               )
          )
          
          # 3. Store only after warm-up and every thin iterations
          if (b > burn_in && (b - burn_in) %% thin == 0) {
               s <- s + 1
               
               BETA_post[s, ] <- beta
               LL_post[s]     <- log_likelihood_probit(beta, y, X)
          }
          
          # 4. Progress
          if (progress && b %% ncat == 0) {
               cat(round(100 * b / B_post, 1), "% completado ...\n", sep = "")
          }
     }
     
     # Return results
     list(
          BETA   = BETA_post,
          LL     = LL_post,
          V_beta = V_beta
     )
}

# Logit model fitting ----------------------------------------------------------

# Hyperparameters
beta0  <- rep(0, p)
Sigma0 <- 10 * diag(1, p)

# Algorithm configuration
B       <- 10000
burn_in <- 10000
thin    <- 10

# Model fitting
muestras_logit <- mcmc_logit(
     y        = y,
     X        = X,
     beta0    = beta0,
     Sigma0   = Sigma0,
     B        = B,
     burn_in  = burn_in,
     thin     = thin,
     seed     = 123,
     progress = TRUE
)

# Acceptance rate
round(muestras_logit$tasa_aceptacion, 3)

# Extract results
BETA_logit <- muestras_logit$BETA
LL_logit   <- muestras_logit$LL
acr_logit  <- muestras_logit$tasa_aceptacion

# Probit model fitting ---------------------------------------------------------

# Hyperparameters
beta0  <- rep(0, p)
Sigma0 <- 10 * diag(1, p)

# Algorithm configuration
B       <- 10000
burn_in <- 10000
thin    <- 10

# Model fitting
muestras_probit <- mcmc_probit(
     y        = y,
     X        = X,
     beta0    = beta0,
     Sigma0   = Sigma0,
     B        = B,
     burn_in  = burn_in,
     thin     = thin,
     seed     = 123,
     progress = TRUE
)

# Extract results
BETA_probit <- muestras_probit$BETA
LL_probit   <- muestras_probit$LL

# Convergence and posterior summary for the logit model ------------------------

# Matrix of posterior quantities: coefficients and log-likelihood
MCMC_logit <- cbind(
     BETA_logit,
     logLik = LL_logit
)

# Effective sample sizes
ESS_logit <- effectiveSize(MCMC_logit)

# Monte Carlo standard errors
MCSE_logit <- apply(
     MCMC_logit,
     2,
     function(x) {
          sd(x) / sqrt(effectiveSize(x))
     }
)

# Posterior summaries
resumen_posterior_logit <- data.frame(
     ess   = ESS_logit,
     mcse  = MCSE_logit,
     media = colMeans(MCMC_logit),
     sd    = apply(MCMC_logit, 2, sd),
     q025  = apply(MCMC_logit, 2, quantile, probs = 0.025),
     q975  = apply(MCMC_logit, 2, quantile, probs = 0.975)
)

round(resumen_posterior_logit, 3)

# Convergence and posterior summary for the probit model -----------------------

# Matrix of posterior quantities: coefficients and log-likelihood
MCMC_probit <- cbind(
     BETA_probit,
     logLik = LL_probit
)

# Effective sample sizes
ESS_probit <- effectiveSize(MCMC_probit)

# Monte Carlo standard errors
MCSE_probit <- apply(
     MCMC_probit,
     2,
     function(x) {
          sd(x) / sqrt(effectiveSize(x))
     }
)

# Posterior summaries
resumen_posterior_probit <- data.frame(
     ess   = ESS_probit,
     mcse  = MCSE_probit,
     media = colMeans(MCMC_probit),
     sd    = apply(MCMC_probit, 2, sd),
     q025  = apply(MCMC_probit, 2, quantile, probs = 0.025),
     q975  = apply(MCMC_probit, 2, quantile, probs = 0.975)
)

round(resumen_posterior_probit, 3)

# Posterior densities of the log-likelihood: logit vs. probit ------------------

pdf(
     file      = paste0("parapolitica_logit_probit_log_verosimilitud.pdf"),
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

dens_LL_logit  <- density(LL_logit)
dens_LL_probit <- density(LL_probit)

ylim_LL <- range(dens_LL_logit$y, dens_LL_probit$y)

plot(
     dens_LL_logit,
     type = "l",
     lwd  = 3,
     col  = 4,
     xlab = "Log-likelihood",
     ylab = "Density",
     ylim = ylim_LL,
     main = ""
)

lines(
     dens_LL_probit,
     lwd = 3,
     col = 2
)

legend(
     "topleft",
     legend = c("Logit", "Probit"),
     col    = c(4, 2),
     lwd    = 3,
     bty    = "n",
     cex    = 1.5
)

dev.off()

# Posterior densities of the coefficients: logit vs. probit -------------------

parametros <- colnames(BETA_logit)

for (j in seq_along(parametros)) {
     pdf(
          file = paste0(
               "parapolitica_logit_probit_posterior_",
               parametros[j],
               ".pdf"
          ),
          width     = 5,
          height    = 5,
          pointsize = 17
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     dens_logit  <- density(BETA_logit[, j])
     dens_probit <- density(BETA_probit[, j])
     
     xlim_j <- range(dens_logit$x, dens_probit$x)
     ylim_j <- range(dens_logit$y, dens_probit$y)
     
     plot(
          dens_logit,
          type = "l",
          lwd  = 3,
          col  = 4,
          xlim = xlim_j,
          ylim = ylim_j,
          main = "",
          xlab = parametros[j],
          ylab = "Density"
     )
     
     lines(
          dens_probit,
          lwd = 3,
          col = 2
     )
     
     box()
     
     dev.off()
}

# Posterior probability curves: logit model -----------------------------------

# Grid of ideal points
ideal_grid <- seq(
     from       = -3,
     to         = 3,
     length.out = 1000
)

# Design matrix for prediction
X_grid <- cbind(
     1,
     ideal = ideal_grid
)

# Posterior probabilities for all iterations
PROB_logit <- sapply(
     seq_len(nrow(BETA_logit)),
     function(b) {
          plogis(X_grid %*% BETA_logit[b, ])
     }
)

# Posterior mean of the curve
PROB_logit_media <- rowMeans(PROB_logit)

# Random selection of 100 posterior curves
set.seed(123)

ind_muestras <- sample(
     x    = seq_len(nrow(BETA_logit)),
     size = 100
)

PROB_logit_muestras <- PROB_logit[, ind_muestras]

# Plot
pdf(
     file      = "parapolitica_logit_curva_probabilidad_muestras.pdf",
     width     = 7,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = ideal_grid,
     y    = PROB_logit_media,
     type = "n",
     ylim = c(0, 1),
     xlab = "Ideal point",
     ylab = "Probability of parapolitics involvement",
     main = ""
)

# Simulated posterior curves
for (j in seq_len(ncol(PROB_logit_muestras))) {
     lines(
          x   = ideal_grid,
          y   = PROB_logit_muestras[, j],
          col = adjustcolor("gray40", 0.18),
          lwd = 1
     )
}

# Posterior mean
lines(
     x   = ideal_grid,
     y   = PROB_logit_media,
     col = 4,
     lwd = 3
)

# Observed data
points(
     x   = parapolitica_modelo$ideal,
     y   = jitter(parapolitica_modelo$parapol, amount = 0.025),
     pch = 16,
     cex = 0.6,
     col = adjustcolor("black", 0.5)
)

box()

dev.off()

# Posterior probability curves: logit vs. probit ------------------------------

# Grid of ideal points
ideal_grid <- seq(
     from       = -3,
     to         = 3,
     length.out = 1000
)

# Design matrix for prediction
X_grid <- cbind(
     1,
     ideal = ideal_grid
)

# Posterior probabilities for all iterations
PROB_logit <- sapply(
     seq_len(nrow(BETA_logit)),
     function(b) {
          plogis(X_grid %*% BETA_logit[b, ])
     }
)

PROB_probit <- sapply(
     seq_len(nrow(BETA_probit)),
     function(b) {
          pnorm(X_grid %*% BETA_probit[b, ])
     }
)

# Posterior summaries of the curves
PROB_logit_media <- rowMeans(PROB_logit)

PROB_logit_IC <- apply(
     PROB_logit,
     1,
     quantile,
     probs = c(0.025, 0.975)
)

PROB_probit_media <- rowMeans(PROB_probit)

PROB_probit_IC <- apply(
     PROB_probit,
     1,
     quantile,
     probs = c(0.025, 0.975)
)

# Plot
pdf(
     file      = "parapolitica_logit_probit_curva_probabilidad.pdf",
     width     = 7,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = ideal_grid,
     y    = PROB_logit_media,
     type = "n",
     ylim = c(0, 1),
     xlab = "Ideal point",
     ylab = "Probability of parapolitics involvement"
)

# 95% posterior bands: logit
lines(
     x   = ideal_grid,
     y   = PROB_logit_IC[1, ],
     col = 4,
     lwd = 1,
     lty = 2
)

lines(
     x   = ideal_grid,
     y   = PROB_logit_IC[2, ],
     col = 4,
     lwd = 1,
     lty = 2
)

# 95% posterior bands: probit
lines(
     x   = ideal_grid,
     y   = PROB_probit_IC[1, ],
     col = 2,
     lwd = 1,
     lty = 2
)

lines(
     x   = ideal_grid,
     y   = PROB_probit_IC[2, ],
     col = 2,
     lwd = 1,
     lty = 2
)

# Posterior means
lines(
     x   = ideal_grid,
     y   = PROB_logit_media,
     col = 4,
     lwd = 3
)

lines(
     x   = ideal_grid,
     y   = PROB_probit_media,
     col = 2,
     lwd = 3
)

# Observed data
points(
     x   = parapolitica_modelo$ideal,
     y   = jitter(parapolitica_modelo$parapol, amount = 0.025),
     pch = 16,
     cex = 0.6,
     col = adjustcolor("black", 0.5)
)

legend(
     "topleft",
     legend = c("Logit", "Probit"),
     col    = c(4, 2),
     lwd    = 3,
     lty    = 1,
     cex    = 1.5,
     bty    = "n"
)

box()

dev.off()

# AUC for each model: logit vs. probit -----------------------------------------

# Function for computing the AUC using the posterior mean of beta
auc_media_posterior <- function(BETA, y, X, link = c("logit", "probit")) {
     require(pROC)
     
     link <- match.arg(link)
     
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     # Posterior mean of the coefficients
     beta_media <- colMeans(BETA)
     
     # Linear predictor
     eta <- as.numeric(X %*% beta_media)
     
     # Estimated probabilities of success
     if (link == "logit") {
          prob_hat <- plogis(eta)
     }
     
     if (link == "probit") {
          prob_hat <- pnorm(eta)
     }
     
     # AUC
     AUC <- as.numeric(
          pROC::auc(
               pROC::roc(
                    response  = y,
                    predictor = prob_hat,
                    levels    = c(0, 1),
                    direction = "<",
                    quiet     = TRUE
               )
          )
     )
     
     list(
          beta_media = beta_media,
          prob_hat   = prob_hat,
          AUC        = AUC
     )
}

AUC_logit_media <- auc_media_posterior(
     BETA = BETA_logit,
     y    = y,
     X    = X,
     link = "logit"
)

AUC_probit_media <- auc_media_posterior(
     BETA = BETA_probit,
     y    = y,
     X    = X,
     link = "probit"
)

# Summary
resumen_AUC_media <- data.frame(
     AUC = c(
          AUC_logit_media$AUC,
          AUC_probit_media$AUC
     )
)

rownames(resumen_AUC_media) <- c("Logit", "Probit")

round(resumen_AUC_media, 3)

# DIC and WAIC for the logit and probit models ---------------------------------

# Stable function for log(1 + exp(x))
log1pexp <- function(x) {
     ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

# Pointwise log-likelihood matrix
# Returns an n x B matrix, where each column corresponds
# to a posterior iteration
log_lik_matrix_binary <- function(BETA, y, X, link = c("logit", "probit")) {
     link <- match.arg(link)
     
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     B <- nrow(BETA)
     n <- length(y)
     
     LL_mat <- matrix(NA, nrow = n, ncol = B)
     
     for (b in seq_len(B)) {
          eta_b <- as.numeric(X %*% BETA[b, ])
          
          if (link == "logit") {
               LL_mat[, b] <- y * eta_b - log1pexp(eta_b)
          }
          
          if (link == "probit") {
               LL_mat[, b] <- ifelse(
                    y == 1,
                    pnorm(eta_b, log.p = TRUE),
                    pnorm(eta_b, lower.tail = FALSE, log.p = TRUE)
               )
          }
     }
     
     LL_mat
}

# Log-likelihood function evaluated at the posterior mean
log_likelihood_binary <- function(beta, y, X, link = c("logit", "probit")) {
     link <- match.arg(link)
     
     y <- as.numeric(y)
     X <- as.matrix(X)
     
     eta <- as.numeric(X %*% beta)
     
     if (link == "logit") {
          ll <- sum(y * eta - log1pexp(eta))
     }
     
     if (link == "probit") {
          ll <- sum(
               ifelse(
                    y == 1,
                    pnorm(eta, log.p = TRUE),
                    pnorm(eta, lower.tail = FALSE, log.p = TRUE)
               )
          )
     }
     
     ll
}

# Auxiliary log-mean-exp function
log_mean_exp <- function(x) {
     m <- max(x)
     m + log(mean(exp(x - m)))
}

# Compute DIC and WAIC
dic_waic_binary <- function(BETA, y, X, link = c("logit", "probit")) {
     link <- match.arg(link)
     
     # Pointwise log-likelihoods
     LL_mat <- log_lik_matrix_binary(
          BETA = BETA,
          y    = y,
          X    = X,
          link = link
     )
     
     # Total log-likelihood by iteration
     LL_total <- colSums(LL_mat)
     
     # DIC
     D_bar <- mean(-2 * LL_total)
     
     beta_bar <- colMeans(BETA)
     
     D_hat <- -2 * log_likelihood_binary(
          beta = beta_bar,
          y    = y,
          X    = X,
          link = link
     )
     
     p_D <- D_bar - D_hat
     DIC <- D_bar + p_D
     
     # WAIC
     lppd <- sum(
          apply(
               LL_mat,
               1,
               log_mean_exp
          )
     )
     
     p_WAIC <- sum(
          apply(
               LL_mat,
               1,
               var
          )
     )
     
     WAIC <- -2 * (lppd - p_WAIC)
     
     data.frame(
          D_bar  = D_bar,
          D_hat  = D_hat,
          p_D    = p_D,
          DIC    = DIC,
          lppd   = lppd,
          p_WAIC = p_WAIC,
          WAIC   = WAIC
     )
}

# Apply to each model
criterios_logit <- dic_waic_binary(
     BETA = BETA_logit,
     y    = y,
     X    = X,
     link = "logit"
)

criterios_probit <- dic_waic_binary(
     BETA = BETA_probit,
     y    = y,
     X    = X,
     link = "probit"
)

# Comparative table
criterios_modelos <- rbind(
     Logit  = criterios_logit,
     Probit = criterios_probit
)

round(criterios_modelos, 3)

# End --------------------------------------------------------------------------