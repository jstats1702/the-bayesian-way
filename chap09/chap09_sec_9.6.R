# Settings ---------------------------------------------------------------------

rm(list = ls())

setwd("~/Dropbox/UN/bayes_book_en")

# Functions --------------------------------------------------------------------

r1 <- function(x) round(x, digits = 1)
r2 <- function(x) round(x, digits = 2)
r3 <- function(x) round(x, digits = 3)
r4 <- function(x) round(x, digits = 4)

stats95 <- function(x) {
     x <- as.matrix(x)
     
     media <- apply(
          X      = x,
          MARGIN = 2,
          FUN    = mean,
          na.rm  = TRUE
     )
     
     sd_post <- apply(
          X      = x,
          MARGIN = 2,
          FUN    = sd,
          na.rm  = TRUE
     )
     
     qnt <- apply(
          X      = x,
          MARGIN = 2,
          FUN    = quantile,
          probs  = c(0.025, 0.975),
          na.rm  = TRUE
     )
     
     stats <- cbind(
          Mean = media,
          SD   = sd_post,
          t(qnt)
     )
     
     colnames(stats) <- c(
          "Mean",
          "SD",
          "2.5%",
          "97.5%"
     )
     
     round(stats, 3)
}

# Data -------------------------------------------------------------------------

# Load the data
sst <- read.table(
     file   = paste0("sst.dat"),
     header = TRUE,
     sep    = ""
)

# Extract geographic location and device type
lon <- sst$lon
lat <- sst$lat
tp0 <- sst$Type

# Extract mean temperature, number of records, and sample size
yi <- sst$temp
ni <- sst$N
Iy <- length(yi)

# Basic descriptive summary of temperature
mean(yi)
sd(yi)

# Recode device type
types <- sort(unique(tp0))
J     <- length(types)

typesc <- c(
     "Bucket",
     "Drifting buoys",
     "ERI",
     "Fixed buoys"
)

# Numeric index associated with device type
tpc <- match(tp0, types)

# Construct indicator variables using Bucket as the reference category
tp <- matrix(
     0,
     nrow = Iy,
     ncol = 3
)

colnames(tp) <- c(
     "d.buoy",
     "eri",
     "f.buoy"
)

for (i in seq_len(3)) {
     tp[tp0 == types[i + 1], i] <- 1
}

# Design matrix: intercept, centered longitude and latitude, and device-type indicators
X <- cbind(
     1,
     lon - mean(lon),
     lat - mean(lat),
     tp
)

p <- ncol(X)

# SST map by intensity and device type -----------------------------------------

nclr <- 7

plotclr <- RColorBrewer::brewer.pal(
     n    = nclr,
     name = "YlOrRd"
)

cls <- classInt::classIntervals(
     var   = sst$temp,
     n     = nclr,
     style = "pretty"
)

colcode <- classInt::findColours(
     clI = cls,
     pal = plotclr
)

# Recode device types
types <- sort(unique(tp0))
J     <- length(types)

typesc <- c(
     "Bucket",
     "Drifting buoys",
     "ERI",
     "Fixed buoys"
)

tpc <- match(tp0, types)

# Symbol codes by device type
pchcode <- c(21, 22, 23, 24)[tpc]

pdf(
     file      = "mediterraneo_eda_tipo_temperatura.pdf",
     width     = 8,
     height    = 8,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3.5, 3.5, 1.5, 1.5),
     mgp   = c(2.0, 0.7, 0),
     bg    = "white"
)

xlim_mapa <- range(lon, na.rm = TRUE) + c(-1.0, 1.0)
ylim_mapa <- range(lat, na.rm = TRUE) + c(-2.5, 2.5)

maps::map(
     database = "world",
     xlim     = xlim_mapa,
     ylim     = ylim_mapa,
     fill     = TRUE,
     col      = "gray85",
     border   = "gray55",
     bg       = "aliceblue",
     xlab     = "Longitude",
     ylab     = "Latitude",
     asp      = 1 / cos(mean(lat, na.rm = TRUE) * pi / 180)
)

box(col = "gray35")

points(
     x   = lon,
     y   = lat,
     pch = pchcode,
     bg  = colcode,
     col = colcode,
     cex = 1.5,
     lwd = 0.5
)

legend(
     "bottomleft",
     legend = names(attr(colcode, "table")),
     fill   = attr(colcode, "palette"),
     border = attr(colcode, "palette"),
     title  = expression("SST (" * degree * C * ")"),
     cex    = 1,
     bty    = "n",
     inset  = 0.02
)

legend(
     "topleft",
     legend = typesc,
     pch    = c(21, 22, 23, 24),
     pt.bg  = "gray90",
     col    = "gray20",
     pt.cex = 1.2,
     cex    = 1.2,
     bty    = "n",
     inset  = 0.02
)

dev.off()

# Linear regression ------------------------------------------------------------

# Center longitude and latitude before fitting
sst$lon_c <- sst$lon - mean(sst$lon, na.rm = TRUE)
sst$lat_c <- sst$lat - mean(sst$lat, na.rm = TRUE)

# Linear model with centered coordinates, device type, and all other
# available covariates, excluding N, Dev, and the original lon and lat
modelo_lm_1 <- lm(
     formula = temp ~ . - N - Dev - lon - lat,
     data    = sst
)

summary(modelo_lm_1)

# Gibbs sampler for the hierarchical SST model ---------------------------------

rigamma <- function(n, shape, scale) {
     1 / rgamma(
          n     = n,
          shape = shape,
          rate  = scale
     )
}

gibbs_tsm <- function(
          y,
          X,
          ni,
          alpha,
          a_sigma,
          b_sigma,
          n_sams,
          n_burn,
          n_skip,
          seed = 123,
          verbose = TRUE
) {
     # Adjustments
     y  <- as.numeric(y)
     X  <- as.matrix(X)
     ni <- as.numeric(ni)
     
     n <- length(y)
     p <- ncol(X)
     
     beta_names <- c(
          "int (bckt)",
          "lon",
          "lat",
          "d.buoy",
          "eri",
          "f.buoy"
     )
     
     colnames(X) <- beta_names
     
     # Total number of iterations
     B_total <- n_burn + n_sams * n_skip
     ncat    <- max(1, floor(0.1 * B_total))
     
     # Fixed quantities
     XtX     <- crossprod(X)
     XtX_inv <- solve(XtX)
     
     # Initial values
     beta <- as.numeric(qr.solve(X, y))
     mu   <- y
     
     resid_mu <- as.numeric(mu - X %*% beta)
     tau2     <- max(var(resid_mu), .Machine$double.eps)
     
     sigma2i <- rep(
          max(var(y), .Machine$double.eps),
          n
     )
     
     sigma2 <- mean(sigma2i)
     
     # Storage
     BETA    <- matrix(NA_real_, nrow = n_sams, ncol = p)
     MU      <- matrix(NA_real_, nrow = n_sams, ncol = n)
     SIGMA2I <- matrix(NA_real_, nrow = n_sams, ncol = n)
     
     TAU    <- numeric(n_sams)
     TAU2   <- numeric(n_sams)
     SIGMA  <- numeric(n_sams)
     SIGMA2 <- numeric(n_sams)
     LL     <- numeric(n_sams)
     
     colnames(BETA)    <- beta_names
     colnames(MU)      <- paste0("mu_", seq_len(n))
     colnames(SIGMA2I) <- paste0("sig2_", seq_len(n))
     
     # Chain
     set.seed(seed)
     
     for (b in seq_len(B_total)) {
          # Update mu_i
          v_mu <- 1 / (
               ni / sigma2i +
                    1 / tau2
          )
          
          m_mu <- v_mu * (
               ni * y / sigma2i +
                    as.numeric(X %*% beta) / tau2
          )
          
          mu <- rnorm(
               n    = n,
               mean = m_mu,
               sd   = sqrt(v_mu)
          )
          
          # Update beta
          V_beta <- tau2 * XtX_inv
          
          m_beta <- as.numeric(
               XtX_inv %*% crossprod(X, mu)
          )
          
          beta <- as.numeric(
               MASS::mvrnorm(
                    n     = 1,
                    mu    = m_beta,
                    Sigma = V_beta
               )
          )
          
          # Update tau2
          resid_mu <- as.numeric(mu - X %*% beta)
          
          a_tau <- (n - 1) / 2
          b_tau <- 0.5 * sum(resid_mu^2)
          
          tau2 <- rigamma(
               n     = 1,
               shape = a_tau,
               scale = b_tau
          )
          
          # Update sigma_i^2
          a_sigma_i <- alpha + 3 / 2
          
          b_sigma_i <- alpha * sigma2 +
               0.5 * ni * (y - mu)^2
          
          sigma2i <- rigamma(
               n     = n,
               shape = a_sigma_i,
               scale = b_sigma_i
          )
          
          # Update sigma^2
          a_sigma_star <- a_sigma + n * (alpha + 1)
          
          b_sigma_star <- b_sigma +
               alpha * sum(1 / sigma2i)
          
          sigma2 <- rgamma(
               n     = 1,
               shape = a_sigma_star,
               rate  = b_sigma_star
          )
          
          # Store posterior sample
          if (b > n_burn && (b - n_burn) %% n_skip == 0) {
               k <- (b - n_burn) / n_skip
               
               BETA[k, ]    <- beta
               MU[k, ]      <- mu
               SIGMA2I[k, ] <- sigma2i
               
               TAU[k]    <- sqrt(tau2)
               TAU2[k]   <- tau2
               SIGMA[k]  <- sqrt(sigma2)
               SIGMA2[k] <- sigma2
               
               LL[k] <- sum(
                    dnorm(
                         x    = y,
                         mean = mu,
                         sd   = sqrt(sigma2i / ni),
                         log  = TRUE
                    )
               )
          }
          
          # Progress
          if (verbose && b %% ncat == 0) {
               cat(
                    sprintf(
                         "%.1f%% completado\n",
                         100 * b / B_total
                    )
               )
          }
     }
     
     # Output
     list(
          BETA    = BETA,
          TAU     = TAU,
          TAU2    = TAU2,
          SIGMA   = SIGMA,
          SIGMA2  = SIGMA2,
          MU      = MU,
          SIGMA2I = SIGMA2I,
          LL      = LL
     )
}

# Model fitting ----------------------------------------------------------------

muestras_tsm <- gibbs_tsm(
     y       = yi,
     X       = X,
     ni      = ni,
     alpha   = 2,
     a_sigma = 1,
     b_sigma = 20,
     n_sams  = 10000,
     n_burn  = 10000,
     n_skip  = 100,
     seed    = 123,
     verbose = TRUE
)

# Number of iterations
B <- length(muestras_tsm$LL)

# Convergence ------------------------------------------------------------------

# Ridge log-likelihood chain
plot(
     muestras_tsm$LL,
     type = "p",
     cex  = 0.3,
     col  = adjustcolor(1, 0.5),
     xlab = "Iteration",
     ylab = "Log-likelihood",
     main = ""
)

dev.off()

# Posterior inference for beta, sigma, and tau ---------------------------------

r3(stats95(muestras_tsm$BETA))
r3(stats95(muestras_tsm$TAU))
r3(stats95(muestras_tsm$SIGMA))

# Posterior intervals for the latent means mu_i --------------------------------

cols <- c(
     "coral2",
     "goldenrod3",
     "forestgreen",
     "slateblue"
)

mu_stat <- stats95(muestras_tsm$MU)[, c("Mean", "2.5%", "97.5%")]

id_dispositivo <- seq_len(nrow(mu_stat))

ylim_mu <- range(
     mu_stat[, c("2.5%", "97.5%")],
     mean(yi, na.rm = TRUE),
     na.rm = TRUE
)

pdf(
     file      = "mediterraneo_intervalos_mu.pdf",
     width     = 7.5,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = id_dispositivo,
     y    = mu_stat[, "Mean"],
     pch  = 20,
     col  = cols[tpc],
     xlab = "Device identifier",
     ylab = expression(mu[i]),
     ylim = ylim_mu,
     lwd  = 1,
     main = ""
)

segments(
     x0     = id_dispositivo,
     y0     = mu_stat[, "2.5%"],
     x1     = id_dispositivo,
     y1     = mu_stat[, "97.5%"],
     col    = cols[tpc],
     length = 0,
     angle  = 90,
     code   = 3
)

abline(
     h   = mean(yi, na.rm = TRUE),
     lty = 2,
     lwd = 2,
     col = "gray30"
)

legend(
     "bottomright",
     legend = typesc,
     pch    = 20,
     col    = cols,
     cex    = 1,
     bty    = "n"
)

dev.off()

# Posterior intervals for the individual variances sigma_i^2 ------------------

cols <- c(
     "coral2",
     "goldenrod3",
     "forestgreen",
     "slateblue"
)

sig2i_stat <- stats95(muestras_tsm$SIGMA2I)[, c("Mean", "2.5%", "97.5%")]

id_dispositivo <- seq_len(nrow(sig2i_stat))

ylim_sig2i <- range(
     sig2i_stat[, c("2.5%", "97.5%")],
     mean(muestras_tsm$SIGMA2, na.rm = TRUE),
     na.rm = TRUE
)

pdf(
     file      = "mediterraneo_intervalos_sigma2.pdf",
     width     = 7.5,
     height    = 5,
     pointsize = 17
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

plot(
     x    = id_dispositivo,
     y    = sig2i_stat[, "Mean"],
     pch  = 20,
     col  = cols[tpc],
     xlab = "Device identifier",
     ylab = expression(sigma[i]^2),
     ylim = ylim_sig2i,
     lwd  = 1,
     main = ""
)

arrows(
     x0     = id_dispositivo,
     y0     = sig2i_stat[, "2.5%"],
     x1     = id_dispositivo,
     y1     = sig2i_stat[, "97.5%"],
     col    = cols[tpc],
     length = 0,
     angle  = 90,
     code   = 3
)

abline(
     h   = mean(muestras_tsm$SIGMA2, na.rm = TRUE),
     lty = 2,
     lwd = 2,
     col = "gray30"
)

dev.off()

# Posterior predictive checks using standardized residuals --------------------

skewness <- function(x) {
     x <- as.numeric(x)
     x <- x[is.finite(x)]
     
     z <- x - mean(x)
     mean(z^3) / sd(x)^3
}

kurtosis_excess <- function(x) {
     x <- as.numeric(x)
     x <- x[is.finite(x)]
     
     z <- x - mean(x)
     mean(z^4) / sd(x)^4 - 3
}

discrepancias <- function(r) {
     c(
          media     = mean(r),
          escala    = sd(r),
          mediana   = median(r),
          iqr       = IQR(r),
          asimetria = skewness(r),
          curtosis  = kurtosis_excess(r)
     )
}

B <- nrow(muestras_tsm$MU)
n <- length(yi)

T_obs <- matrix(NA_real_, nrow = B, ncol = 6)
T_rep <- matrix(NA_real_, nrow = B, ncol = 6)

colnames(T_obs) <- colnames(T_rep) <- c(
     "media",
     "escala",
     "mediana",
     "iqr",
     "asimetria",
     "curtosis"
)

set.seed(123)

for (b in seq_len(B)) {
     mu_b  <- muestras_tsm$MU[b, ]
     sig_b <- sqrt(muestras_tsm$SIGMA2I[b, ])
     
     y_rep <- rnorm(
          n    = n,
          mean = mu_b,
          sd   = sig_b / sqrt(ni)
     )
     
     r_obs <- (yi - mu_b) / (sig_b / sqrt(ni))
     r_rep <- (y_rep - mu_b) / (sig_b / sqrt(ni))
     
     T_obs[b, ] <- discrepancias(r_obs)
     T_rep[b, ] <- discrepancias(r_rep)
}

ppp <- colMeans(T_rep <= T_obs)

round(ppp, 3)

# Histograms of the posterior predictive discrepancies ------------------------

n_estadisticos <- ncol(T_rep)

nombres_estadisticos <- colnames(T_rep)

etiquetas_estadisticos <- c(
     media     = "Overall centering",
     escala    = "Scale",
     mediana   = "Median",
     iqr       = "Interquartile range",
     asimetria = "Skewness",
     curtosis  = "Kurtosis"
)

for (j in seq_len(n_estadisticos)) {
     nombre_j <- nombres_estadisticos[j]
     
     etiqueta_j <- ifelse(
          nombre_j %in% names(etiquetas_estadisticos),
          etiquetas_estadisticos[nombre_j],
          nombre_j
     )
     
     x_j <- T_rep[, j]
     x_j <- x_j[is.finite(x_j)]
     
     archivo_j <- paste0(
          "mediterraneo_histograma_ppp_",
          nombre_j,
          ".pdf"
     )
     
     pdf(
          file      = archivo_j,
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
          x      = x_j,
          breaks = 25,
          freq   = FALSE,
          col    = "gray85",
          border = "white",
          xlab   = expression(t^{rep}),
          ylab   = "Density",
          main   = ""
     )
     
     t_obs_j <- mean(
          T_obs[, j],
          na.rm = TRUE
     )
     
     abline(
          v   = t_obs_j,
          lty = 2,
          lwd = 2,
          col = "gray20"
     )
     
     legend(
          "topright",
          legend = paste0(
               "ppp = ",
               round(ppp[j], 3)
          ),
          bty = "n",
          cex = 1.2
     )
     
     dev.off()
}

# Basic residual check ---------------------------------------------------------

# Posterior means of the device-specific parameters
mu_hat <- colMeans(
     muestras_tsm$MU,
     na.rm = TRUE
)

sig2i_hat <- colMeans(
     muestras_tsm$SIGMA2I,
     na.rm = TRUE
)

# Raw and standardized residuals
resid_crudo <- yi - mu_hat

resid_est <- resid_crudo / sd(resid_crudo)

# Colors by device type
cols <- c(
     "coral2",
     "goldenrod3",
     "forestgreen",
     "slateblue"
)

# Histogram of standardized residuals
pdf(
     file      = "mediterraneo_residuales_histograma.pdf",
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
     x           = resid_est,
     breaks      = 20,
     probability = TRUE,
     col         = "gray85",
     border      = "white",
     xlim        = 4 * c(-1, 1),
     xlab        = "Standardized residual",
     ylab        = "Density",
     main        = ""
)

curve(
     dnorm(x),
     add = TRUE,
     lwd = 2,
     col = "gray30"
)

abline(
     v   = 0,
     lty = 2,
     lwd = 2,
     col = "gray30"
)

box()

dev.off()

# Normal Q-Q plot
pdf(
     file      = "mediterraneo_residuales_qqplot.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

qqnorm(
     resid_est,
     pch  = 16,
     cex  = 1.5,
     col  = adjustcolor("gray30", 0.7),
     xlim = 4 * c(-1, 1),
     ylim = 4 * c(-1, 1),
     xlab = "Theoretical quantiles",
     ylab = "Sample quantiles",
     main = ""
)

qqline(
     resid_est,
     lwd = 2,
     lty = 2,
     col = "gray30"
)

dev.off()

# Residuals versus fitted values
pdf(
     file      = "mediterraneo_residuales_ajustados.pdf",
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
     x    = mu_hat,
     y    = resid_est,
     pch  = 16,
     cex  = 1.1,
     col  = adjustcolor(cols[tpc], 0.7),
     xlab = expression(hat(mu)[i]),
     ylab = "Standardized residual",
     main = ""
)

abline(
     h   = 0,
     lty = 2,
     lwd = 2,
     col = "gray30"
)

legend(
     "bottomright",
     legend = typesc,
     pch    = 20,
     col    = cols,
     cex    = 1,
     bty    = "n"
)

dev.off()

# Residuals by device identifier
pdf(
     file      = "mediterraneo_residuales_dispositivo.pdf",
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
     x    = seq_along(resid_est),
     y    = resid_est,
     pch  = 16,
     cex  = 1.1,
     col  = adjustcolor(cols[tpc], 0.7),
     xlab = "Device identifier",
     ylab = "Standardized residual",
     main = ""
)

abline(
     h   = 0,
     lty = 2,
     lwd = 2,
     col = "gray30"
)

legend(
     "bottomright",
     legend = typesc,
     pch    = 20,
     col    = cols,
     cex    = 1,
     bty    = "n"
)

dev.off()

# Posterior prediction over a geographic grid ---------------------------------

# Grid sizes
nsamp1 <- 100
nsamp2 <- 100

lon_grid <- seq(
     from       = min(lon, na.rm = TRUE),
     to         = max(lon, na.rm = TRUE),
     length.out = nsamp1
)

lat_grid <- seq(
     from       = min(lat, na.rm = TRUE),
     to         = max(lat, na.rm = TRUE),
     length.out = nsamp2
)

lon_grid_c <- lon_grid - mean(lon, na.rm = TRUE)
lat_grid_c <- lat_grid - mean(lat, na.rm = TRUE)

grid_pred <- expand.grid(
     lon_c = lon_grid_c,
     lat_c = lat_grid_c
)

predecir_tsm_grilla <- function(
          tp_new,
          muestras,
          alpha,
          n_new,
          seed = 123
) {
     B      <- nrow(muestras$BETA)
     n_grid <- nrow(grid_pred)
     
     X_new <- cbind(
          1,
          grid_pred$lon_c,
          grid_pred$lat_c,
          matrix(
               tp_new,
               nrow  = n_grid,
               ncol  = length(tp_new),
               byrow = TRUE
          )
     )
     
     colnames(X_new) <- colnames(muestras$BETA)
     
     suma_y  <- rep(0, n_grid)
     suma_y2 <- rep(0, n_grid)
     
     for (b in seq_len(B)) {
          media_reg_b <- as.numeric(
               X_new %*% muestras$BETA[b, ]
          )
          
          mu_new_b <- rnorm(
               n    = n_grid,
               mean = media_reg_b,
               sd   = muestras$TAU[b]
          )
          
          sig2_new_b <- rigamma(
               n     = n_grid,
               shape = alpha + 1,
               scale = alpha * muestras$SIGMA2[b]
          )
          
          y_new_b <- rnorm(
               n    = n_grid,
               mean = mu_new_b,
               sd   = sqrt(sig2_new_b / n_new)
          )
          
          suma_y  <- suma_y + y_new_b
          suma_y2 <- suma_y2 + y_new_b^2
     }
     
     media_pred <- suma_y / B
     
     var_pred <- suma_y2 / B - media_pred^2
     var_pred <- pmax(var_pred, 0)
     
     list(
          media = matrix(
               media_pred,
               nrow = nsamp1,
               ncol = nsamp2
          ),
          sd = matrix(
               sqrt(var_pred),
               nrow = nsamp1,
               ncol = nsamp2
          )
     )
}

# Prediction for each device type
tipos_pred <- list(
     bucket = c(0, 0, 0),
     d.buoy = c(1, 0, 0),
     eri    = c(0, 1, 0),
     f.buoy = c(0, 0, 1)
)

etiquetas_tipos_pred <- c(
     bucket = "Bucket",
     d.buoy = "Drifting buoys",
     eri    = "ERI",
     f.buoy = "Fixed buoys"
)

pred_tsm <- lapply(
     names(tipos_pred),
     function(nombre_tipo) {
          predecir_tsm_grilla(
               tp_new   = tipos_pred[[nombre_tipo]],
               muestras = muestras_tsm,
               alpha    = 2,
               n_new    = 1,
               seed     = 123
          )
     }
)

names(pred_tsm) <- names(tipos_pred)

# Function for plotting predictive maps by device type
graficar_prediccion_tsm <- function(
          pred,
          variable = c("media", "sd")
) {
     zlim <- range(
          unlist(
               lapply(
                    pred,
                    function(x) {
                         as.vector(x[[variable]])
                    }
               )
          ),
          na.rm = TRUE
     )
     
     archivos <- character(length(pred))
     names(archivos) <- names(pred)
     
     for (nombre_tipo in names(pred)) {
          z_j <- pred[[nombre_tipo]][[variable]]
          
          nombre_archivo <- gsub(
               pattern     = "[^A-Za-z0-9]+",
               replacement = "_",
               x           = nombre_tipo
          )
          
          archivo <- paste0(
               "mediterraneo_prediccion",
               "_",
               variable,
               "_tsm_",
               nombre_archivo,
               ".pdf"
          )
          
          archivos[nombre_tipo] <- archivo
          
          pdf(
               file      = archivo,
               width     = 6,
               height    = 5,
               pointsize = 15
          )
          
          oldpar <- par(no.readonly = TRUE)
          on.exit(par(oldpar), add = TRUE)
          
          par(
               mfrow = c(1, 1),
               mar   = c(3.0, 3.0, 1.4, 1.4),
               mgp   = c(1.75, 0.75, 0)
          )
          
          fields::image.plot(
               x          = lon_grid,
               y          = lat_grid,
               z          = z_j,
               zlim       = zlim,
               col        = fields::tim.colors(64),
               xlab       = "Longitude",
               ylab       = "Latitude",
               main       = "",
               legend.lab = "",
               asp        = 1 / cos(mean(lat, na.rm = TRUE) * pi / 180)
          )
          
          xlim_mapa <- range(lon, na.rm = TRUE) + c(-1.0, 1.0)
          ylim_mapa <- range(lat, na.rm = TRUE) + c(-2.5, 2.5)
          
          mapa <- maps::map(
               database = "world",
               xlim     = xlim_mapa,
               ylim     = ylim_mapa,
               asp      = 1 / cos(mean(lat, na.rm = TRUE) * pi / 180),
               plot     = F
          )
          
          lines(
               x   = mapa$x,
               y   = mapa$y,
               col = "black",
               lwd = 2
          )
          
          points(
               x   = lon,
               y   = lat,
               pch = 20,
               cex = 1,
               col = "gray20"
          )
          
          box()
          
          dev.off()
     }
     
     invisible(archivos)
}

# Posterior predictive mean maps
archivos_media <- graficar_prediccion_tsm(
     pred     = pred_tsm,
     variable = "media"
)

# End --------------------------------------------------------------------------