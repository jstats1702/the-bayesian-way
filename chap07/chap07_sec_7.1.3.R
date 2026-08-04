# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
setwd("~/Dropbox/UN/bayes_book_en")

# Packages
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(sf)))
suppressMessages(suppressWarnings(library(stringi)))

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
dat <- dat[idx_muestra, , drop = FALSE]

dim(dat)
table(dat$estu_depto_reside)

# Sufficient statistics --------------------------------------------------------

# m : number of groups (departments)
deptos <- sort(unique(dat$estu_cod_reside_depto))
m      <- length(deptos)

# n : number of individuals (students)
n <- nrow(dat)

# Data processing
# y  : student scores (c)
# Y  : student scores by department (list)
# g  : sequential department identifier (c)
# nj : number of students by department (c)
# yb : means by department (c)
# s2 : variances by department (c)

y <- dat$punt_matematicas
Y <- vector(mode = "list", length = m)
g <- rep(NA_integer_, n)

for (j in seq_len(m)) {
     idx    <- dat$estu_cod_reside_depto == deptos[j]
     g[idx] <- j
     Y[[j]] <- y[idx]
}

# Table
estadisticos <- dat %>%
     group_by(estu_cod_reside_depto) %>%
     summarise(
          codigo  = first(estu_cod_reside_depto),
          nombre  = first(estu_depto_reside),
          nj      = dplyr::n(),
          yb      = mean(punt_matematicas, na.rm = TRUE),
          med     = median(punt_matematicas, na.rm = TRUE),
          s2      = var(punt_matematicas, na.rm = TRUE),
          s       = sd(punt_matematicas, na.rm = TRUE),
          min     = min(punt_matematicas, na.rm = TRUE),
          max     = max(punt_matematicas, na.rm = TRUE),
          .groups = "drop"
     ) %>%
     arrange(codigo) %>%
     select(codigo, nombre, nj, yb, med, s2, s, min, max)

as.data.frame(estadisticos[, c(2, 3, 4, 5, 7, 8, 9)])

# Summary vectors
nj <- estadisticos$nj
yb <- estadisticos$yb
s2 <- estadisticos$s2

# Overall mean
round(mean(y), 3)

# Map with average scores ------------------------------------------------------

# Read GeoJSON of Colombian departments
url_geojson <- paste0(
     "https://raw.githubusercontent.com/caticoa3/colombia_mapa/master/",
     "co_2018_MGN_DPTO_POLITICO.geojson"
)

map_col <- sf::st_read(url_geojson, quiet = TRUE)

# Prepare map names
map_col <- map_col %>%
     mutate(
          nombre_key = DPTO_CNMBR %>%
               stringi::stri_trans_general("Latin-ASCII") %>%
               toupper()
     )

# Prepare departmental data
estadisticos_mapa <- estadisticos %>%
     mutate(
          nombre_key = nombre %>%
               stringi::stri_trans_general("Latin-ASCII") %>%
               toupper(),
          nombre_key = case_when(
               nombre_key == "BOGOTA"          ~ "BOGOTA, D.C.",
               nombre_key == "VALLE"           ~ "VALLE DEL CAUCA",
               nombre_key == "NORTE SANTANDER" ~ "NORTE DE SANTANDER",
               TRUE                            ~ nombre_key
          )
     )

# Join map and data
map_col_yb <- map_col %>%
     left_join(
          estadisticos_mapa,
          by = "nombre_key"
     )

# Review unmatched departments
map_col_yb %>%
     filter(is.na(yb)) %>%
     select(DPTO_CNMBR, nombre_key)

# Exclude San Andrés and Providencia to avoid empty space on the left
map_col_yb_cont <- map_col_yb %>%
     filter(
          !grepl(
               pattern     = "SAN ANDRES|PROVIDENCIA|SANTA CATALINA",
               x           = nombre_key,
               ignore.case = TRUE
          )
     )

# Map limits
bbox <- sf::st_bbox(map_col_yb_cont)

pad_x <- 0.35
pad_y <- 0.35

# Map
p_mapa <- ggplot(map_col_yb_cont) +
     geom_sf(
          aes(fill = yb),
          color     = "gray40",
          linewidth = 0.15
     ) +
     scale_fill_gradientn(
          colors   = c("#f3e8ff", "#d8b4fe", "#c084fc", "#a855f7", "#7e22ce"),
          na.value = "gray90",
          name     = expression(bar(y)[j])
     ) +
     coord_sf(
          xlim   = c(bbox["xmin"] - pad_x, bbox["xmax"] + pad_x),
          ylim   = c(bbox["ymin"] - pad_y, bbox["ymax"] + pad_y),
          expand = FALSE
     ) +
     labs(
          x = "Longitude",
          y = "Latitude"
     ) +
     theme_bw(base_size = 13) +
     theme(
          legend.position = "right",
          panel.border = element_rect(
               color     = "black",
               fill      = NA,
               linewidth = 0.6
          ),
          panel.grid.major = element_line(
               color     = "gray85",
               linewidth = 0.25
          ),
          panel.grid.minor = element_blank(),
          axis.title       = element_text(size = 12),
          axis.text        = element_text(size = 10),
          plot.title       = element_text(hjust = 0),
          plot.subtitle    = element_text(hjust = 0),
          plot.caption     = element_text(hjust = 1),
          plot.margin      = margin(t = 6, r = 6, b = 6, l = 6)
     )

p_mapa

# Save the map
ggsave(
     filename = "matematicas_modelo_jerarquico_normal_mapa_promedio.pdf",
     plot     = p_mapa,
     width    = 5.8,
     height   = 7.2,
     units    = "in"
)

dev.off()

# Frequentist ranking ----------------------------------------------------------

pdf(
     file      = "matematicas_modelo_jerarquico_normal_ranking_promedio.pdf",
     width     = 5.8,
     height    = 7.2,
     pointsize = 13
)

# Ranking based on the sample mean
par(
     mfrow = c(1, 1),
     mar   = c(4, 9.5, 1.5, 1),
     mgp   = c(2.5, 0.75, 0)
)

ord <- order(yb)

# Violet tones
col_puntos <- adjustcolor("#7e22ce", alpha.f = 0.5)
col_media  <- "#7e22ce"
col_linea  <- adjustcolor("#a855f7", alpha.f = 0.45)
col_ref    <- adjustcolor("gray", alpha.f = 0.85)

# Initialize empty plot
plot(
     x    = c(0, 100),
     y    = c(0.5, m + 0.5),
     type = "n",
     xlab = "Score",
     ylab = "",
     main = "",
     yaxt = "n",
     yaxs = "i"
)

# Horizontal reference lines
abline(
     h   = 1:m,
     col = "lightgray",
     lwd = 1
)

# Vertical line at the national mean
abline(
     v   = 50,
     col = col_ref,
     lwd = 3
)

# Individual points by department according to the mean ranking
for (l in 1:m) {
     j <- ord[l]
     
     points(
          x   = Y[[j]],
          y   = rep(l, nj[j]),
          pch = 16,
          cex = 0.4,
          col = col_puntos
     )
}

# Connect sample means
lines(
     x    = yb[ord],
     y    = 1:m,
     type = "l",
     col  = col_linea,
     lwd  = 1.5
)

# Add sample mean points
points(
     x   = yb[ord],
     y   = 1:m,
     pch = 16,
     cex = 1.1,
     col = col_media
)

# Labels with department names
axis(
     side   = 2,
     at     = 1:m,
     labels = estadisticos$nombre[ord],
     las    = 2
)

dev.off()

# Histogram of the means -------------------------------------------------------

pdf(
     file      = "matematicas_modelo_jerarquico_normal_histograma_promedios.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Histogram of the group means
hist(
     x      = yb,
     freq   = FALSE,
     breaks = 15,
     col    = adjustcolor("#d8b4fe", alpha.f = 0.65),
     border = "white",
     xlab   = "Mean",
     ylab   = "Density",
     main   = ""
)

abline(
     v   = mean(y),
     col = "gray",
     lwd = 3
)

dev.off()

# Scatter plot of the mean versus size ----------------------------------------

pdf(
     file      = "matematicas_modelo_jerarquico_normal_promedio_vs_tamano.pdf",
     width     = 5,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Scatter plot: group size vs. mean
plot(
     x    = nj,
     y    = yb,
     xlab = "Group size",
     ylab = "Mean",
     pch  = 16,
     cex  = 1.2,
     col  = adjustcolor("#7e22ce", alpha.f = 0.60)
)

abline(
     h   = mean(y, na.rm = TRUE),
     col = "gray",
     lwd = 3
)

dev.off()

# Prior distribution ----------------------------------------------------------

# Hyperparameters
mu0 <- 50
g20 <- 10^2

eta0 <- 1
t20  <- 10^2

nu0 <- 1
s20 <- 10^2

# Gibbs sampler ----------------------------------------------------------------

mcmc <- function(B, y, nj, yb, s2, mu0, g20, eta0, t20, nu0, s20) {
     # Frequency for printing progress
     ncat <- max(1, floor(B / 10))
     
     # Total number of observations and groups
     n <- sum(nj)
     m <- length(nj)
     
     # Chain storage
     THETA <- matrix(NA_real_, nrow = B, ncol = m + 4)
     
     # Initial values
     theta <- yb
     sig2  <- mean(s2, na.rm = TRUE)
     mu    <- mean(theta, na.rm = TRUE)
     tau2  <- var(theta, na.rm = TRUE)
     
     # MCMC chain
     for (b in seq_len(B)) {
          # Update theta_j
          v_theta <- 1 / (1 / tau2 + nj / sig2)
          m_theta <- v_theta * (mu / tau2 + nj * yb / sig2)
          theta   <- rnorm(n = m, mean = m_theta, sd = sqrt(v_theta))
          
          # Update sigma^2
          ss_sig2 <- sum((nj - 1) * s2 + nj * (yb - theta)^2)
          a_sig2  <- 0.5 * (nu0 + n)
          b_sig2  <- 0.5 * (nu0 * s20 + ss_sig2)
          sig2    <- 1 / rgamma(n = 1, shape = a_sig2, rate = b_sig2)
          
          # Update mu
          v_mu <- 1 / (1 / g20 + m / tau2)
          m_mu <- v_mu * (mu0 / g20 + sum(theta) / tau2)
          mu   <- rnorm(n = 1, mean = m_mu, sd = sqrt(v_mu))
          
          # Update tau^2
          ss_tau2 <- sum((theta - mu)^2)
          a_tau2  <- 0.5 * (eta0 + m)
          b_tau2  <- 0.5 * (eta0 * t20 + ss_tau2)
          tau2    <- 1 / rgamma(n = 1, shape = a_tau2, rate = b_tau2)
          
          # Log-likelihood using sufficient statistics by group
          ll <- sum(
               -0.5 * nj * log(2 * pi * sig2) -
                    0.5 * ((nj - 1) * s2 + nj * (yb - theta)^2) / sig2
          )
          
          # Store iteration
          THETA[b, ] <- c(theta, sig2, mu, tau2, ll)
          
          # Print progress
          if (b %% ncat == 0) {
               cat(100 * round(b / B, 1), "% completado ... \n", sep = "")
          }
     }
     
     # Output
     colnames(THETA) <- c(paste0("theta", seq_len(m)), "sig2", "mu", "tau2", "ll")
     THETA <- as.data.frame(THETA)
     
     return(list(THETA = THETA))
}

# Model fitting ----------------------------------------------------------------

set.seed(123)
chain <- mcmc(B = 10000, y, nj, yb, s2, mu0, g20, eta0, t20, nu0, s20)

# Function for plotting chains -------------------------------------------------

plot_cadena <- function(x, ylab, col, cex, file) {
     pdf(
          file      = file,
          width     = 6,
          height    = 4,
          pointsize = 15
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

# Convergence ------------------------------------------------------------------

# Chains for the log-likelihood
plot_cadena(
     x    = chain$THETA$ll,
     ylab = "Log-likelihood",
     col  = 1,
     cex  = 0.5,
     file = "matematicas_modelo_jerarquico_normal_cadena_logverosimilitud.pdf"
)

# Effective sample sizes
neff <- coda::effectiveSize(
     coda::as.mcmc(chain$THETA)
)

# Monte Carlo standard error
EEMC <- apply(
     X      = chain$THETA,
     MARGIN = 2,
     FUN    = sd,
     na.rm  = TRUE
) / sqrt(neff)

# Monte Carlo coefficient of variation
media_post <- colMeans(
     chain$THETA,
     na.rm = TRUE
)

CVMC <- EEMC / abs(media_post)

round(summary(neff), digits = 1)
round(summary(EEMC), digits = 4)
round(summary(CVMC), digits = 4)

# Plot posterior ---------------------------------------------------------------

plot_posterior <- function(x, xlab, file, legend = TRUE, legend_pos = "topleft") {
     pdf(
          file      = file,
          width     = 5,
          height    = 5,
          pointsize = 15
     )
     
     on.exit(dev.off())
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     hist(
          x      = x,
          freq   = FALSE,
          breaks = 30,
          col    = adjustcolor("gray90", alpha.f = 0.65),
          border = "white",
          xlab   = xlab,
          ylab   = "Density",
          main   = ""
     )
     
     abline(
          v   = mean(x),
          col = 2,
          lwd = 2,
          lty = 2
     )
     
     abline(
          v   = quantile(x, probs = c(0.025, 0.975)),
          col = 4,
          lwd = 2,
          lty = 2
     )
     
     if (legend) {
          legend(
               x      = legend_pos,
               legend = c("Mean", "95% CrI"),
               col    = c(2, 4),
               fill   = c(2, 4),
               border = c(2, 4),
               bty    = "n"
          )
     }
}

# Inference --------------------------------------------------------------------

# Chains for eta, mu, sigma, and tau
PAR <- cbind(
     eta   = chain$THETA$sig2 / (chain$THETA$sig2 + chain$THETA$tau2),
     mu    = chain$THETA$mu,
     sigma = sqrt(chain$THETA$sig2),
     tau   = sqrt(chain$THETA$tau2)
)

# Posterior mu
plot_posterior(
     x          = PAR[, "mu"],
     xlab       = expression(mu),
     file       = "matematicas_modelo_jerarquico_normal_posterior_mu.pdf",
     legend     = TRUE,
     legend_pos = "topleft"
)

# Posterior tau
plot_posterior(
     x          = PAR[, "tau"],
     xlab       = expression(tau),
     file       = "matematicas_modelo_jerarquico_normal_posterior_tau.pdf",
     legend     = FALSE,
     legend_pos = "topleft"
)

# Posterior sigma
plot_posterior(
     x          = PAR[, "sigma"],
     xlab       = expression(sigma),
     file       = "matematicas_modelo_jerarquico_normal_posterior_sigma.pdf",
     legend     = FALSE,
     legend_pos = "topleft"
)

# Posterior eta
plot_posterior(
     x          = PAR[, "eta"],
     xlab       = expression(eta),
     file       = "matematicas_modelo_jerarquico_normal_posterior_eta.pdf",
     legend     = FALSE,
     legend_pos = "topleft"
)

# Posterior summary: mean, CV (%), and quantiles
tab <- cbind(
     `Media`  = colMeans(PAR),
     `CV`     = abs(apply(PAR, 2, sd) / colMeans(PAR)),
     `Q2.5%`  = apply(PAR, 2, quantile, probs = 0.025),
     `Q97.5%` = apply(PAR, 2, quantile, probs = 0.975)
)

round(tab, 3)

# Bayesian ranking -------------------------------------------------------------

pdf(
     file      = "matematicas_modelo_jerarquico_normal_ranking_bayesiano.pdf",
     width     = 5.8,
     height    = 7.2,
     pointsize = 13
)

# Bayesian ranking
par(
     mfrow = c(1, 1),
     mar   = c(4, 9.5, 1.5, 1),
     mgp   = c(2.5, 0.75, 0)
)

# Posterior summaries
THETA_theta <- as.matrix(chain$THETA[, seq_len(m)])

ids2 <- estadisticos$nombre
that <- colMeans(THETA_theta)

ic1 <- apply(
     X      = THETA_theta,
     MARGIN = 2,
     FUN    = quantile,
     probs  = c(0.025, 0.975)
)

# Sort by posterior mean
ord <- order(that)

ids2 <- ids2[ord]
that <- that[ord]
ic1  <- ic1[, ord]

# Colors according to the interval position relative to 50
colo <- rep(2, m)
colo[ic1[1, ] > 50] <- 1
colo[ic1[2, ] < 50] <- 3

colo <- c("royalblue", "black", "red")[colo]

# Reference color
col_ref <- adjustcolor("gray", alpha.f = 0.85)

# Initialize empty plot
plot(
     x    = c(0, 100),
     y    = c(0.5, m + 0.5),
     type = "n",
     xlab = "Score",
     ylab = "",
     main = "",
     yaxt = "n",
     yaxs = "i"
)

# Horizontal reference lines
abline(
     h   = 1:m,
     col = "lightgray",
     lwd = 1
)

# Vertical line at the reference mean
abline(
     v   = 50,
     col = col_ref,
     lwd = 3
)

# Credible intervals and posterior means
for (j in 1:m) {
     segments(
          x0  = ic1[1, j],
          y0  = j,
          x1  = ic1[2, j],
          y1  = j,
          col = colo[j],
          lwd = 1.5
     )
     
     points(
          x   = that[j],
          y   = j,
          pch = 16,
          cex = 1.1,
          col = colo[j]
     )
}

# Labels with department names
axis(
     side   = 2,
     at     = 1:m,
     labels = ids2,
     las    = 2
)

dev.off()

# CV of the theta values -------------------------------------------------------

# Posterior coefficient of variation of theta_j
that <- apply(
     X      = chain$THETA[, seq_len(m)],
     MARGIN = 2,
     FUN    = mean,
     na.rm  = TRUE
)

shat <- apply(
     X      = chain$THETA[, seq_len(m)],
     MARGIN = 2,
     FUN    = sd,
     na.rm  = TRUE
)

cv_b <- abs(shat / that)

round(summary(cv_b), 3)

# Map with posterior means of theta --------------------------------------------

# Read GeoJSON of Colombian departments
url_geojson <- paste0(
     "https://raw.githubusercontent.com/caticoa3/colombia_mapa/master/",
     "co_2018_MGN_DPTO_POLITICO.geojson"
)

map_col <- sf::st_read(url_geojson, quiet = TRUE)

# Posterior mean of each theta_j
theta_cols <- paste0("theta", seq_len(m))

theta_post <- colMeans(
     as.matrix(chain$THETA[, theta_cols]),
     na.rm = TRUE
)

# Prepare map names
map_col <- map_col %>%
     mutate(
          nombre_key = DPTO_CNMBR %>%
               stringi::stri_trans_general("Latin-ASCII") %>%
               toupper()
     )

# Prepare departmental data
estadisticos_mapa <- estadisticos %>%
     mutate(
          theta_post = as.numeric(theta_post),
          nombre_key = nombre %>%
               stringi::stri_trans_general("Latin-ASCII") %>%
               toupper(),
          nombre_key = case_when(
               nombre_key == "BOGOTA"          ~ "BOGOTA, D.C.",
               nombre_key == "VALLE"           ~ "VALLE DEL CAUCA",
               nombre_key == "NORTE SANTANDER" ~ "NORTE DE SANTANDER",
               TRUE                            ~ nombre_key
          )
     )

# Join map and data
map_col_theta <- map_col %>%
     left_join(
          estadisticos_mapa,
          by = "nombre_key"
     )

# Review unmatched departments
map_col_theta %>%
     filter(is.na(theta_post)) %>%
     select(DPTO_CNMBR, nombre_key)

# Exclude San Andrés and Providencia to avoid empty space on the left
map_col_theta_cont <- map_col_theta %>%
     filter(
          !grepl(
               pattern     = "SAN ANDRES|PROVIDENCIA|SANTA CATALINA",
               x           = nombre_key,
               ignore.case = TRUE
          )
     )

# Map limits
bbox <- sf::st_bbox(map_col_theta_cont)

pad_x <- 0.35
pad_y <- 0.35

# Map
p_mapa <- ggplot(map_col_theta_cont) +
     geom_sf(
          aes(fill = theta_post),
          color     = "gray40",
          linewidth = 0.15
     ) +
     scale_fill_gradientn(
          colors   = c("#f3e8ff", "#d8b4fe", "#c084fc", "#a855f7", "#7e22ce"),
          na.value = "gray90",
          name     = expression(hat(theta)[j])
     ) +
     coord_sf(
          xlim   = c(bbox["xmin"] - pad_x, bbox["xmax"] + pad_x),
          ylim   = c(bbox["ymin"] - pad_y, bbox["ymax"] + pad_y),
          expand = FALSE
     ) +
     labs(
          x = "Longitude",
          y = "Latitude"
     ) +
     theme_bw(base_size = 13) +
     theme(
          legend.position = "right",
          panel.border = element_rect(
               color     = "black",
               fill      = NA,
               linewidth = 0.6
          ),
          panel.grid.major = element_line(
               color     = "gray85",
               linewidth = 0.25
          ),
          panel.grid.minor = element_blank(),
          axis.title       = element_text(size = 12),
          axis.text        = element_text(size = 10),
          plot.title       = element_text(hjust = 0),
          plot.subtitle    = element_text(hjust = 0),
          plot.caption     = element_text(hjust = 1),
          plot.margin      = margin(t = 6, r = 6, b = 6, l = 6)
     )

p_mapa

# Save the map
ggsave(
     filename = "matematicas_modelo_jerarquico_normal_mapa_media_posterior_theta.pdf",
     plot     = p_mapa,
     width    = 5.8,
     height   = 7.2,
     units    = "in"
)

dev.off()

# Shrinkage --------------------------------------------------------------------

# Posterior estimate vs. sample mean
pdf(
     file      = "matematicas_modelo_jerarquico_normal_theta_hat_vs_promedio.pdf",
     width     = 5,
     height    = 5,
     pointsize = 20
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3.3, 1.4, 1.4),
     mgp   = c(2, 0.75, 0)
)

# Posterior means of theta_j
theta_hat <- colMeans(
     as.matrix(chain$THETA[, seq_len(m)]),
     na.rm = TRUE
)

plot(
     x    = yb,
     y    = theta_hat,
     xlim = range(yb, theta_hat, na.rm = TRUE),
     ylim = range(yb, theta_hat, na.rm = TRUE),
     xlab = expression(bar(italic(y))[j]),
     ylab = expression(hat(theta)[j]),
     main = "",
     pch  = 16,
     cex  = 1.2,
     col  = adjustcolor("#7e22ce", alpha.f = 0.60)
)

abline(
     a   = 0,
     b   = 1,
     col = "gray",
     lwd = 3
)

dev.off()

# Difference between sample mean and posterior mean vs. sample size
pdf(
     file      = "matematicas_modelo_jerarquico_normal_theta_hat_vs_tamano.pdf",
     width     = 5,
     height    = 5,
     pointsize = 20
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3.3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

d_theta <- yb - theta_hat

plot(
     x    = nj,
     y    = d_theta,
     xlim = range(nj, na.rm = TRUE),
     ylim = c(-1, 1) * max(abs(d_theta), na.rm = TRUE),
     xlab = "Group size",
     ylab = expression(bar(italic(y))[j] - hat(theta)[j]),
     main = "",
     pch  = 16,
     cex  = 1.2,
     col  = adjustcolor("#7e22ce", alpha.f = 0.60)
)

abline(
     h   = 0,
     col = "gray",
     lwd = 3
)

dev.off()

# Shrinkage: visual comparison between theta_hat and the sample mean
pdf(
     file      = "matematicas_modelo_jerarquico_normal_contraccion_theta_promedio.pdf",
     width     = 5,
     height    = 5,
     pointsize = 20
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Posterior means of theta_j
theta_hat <- colMeans(
     as.matrix(chain$THETA[, seq_len(m)]),
     na.rm = TRUE
)

# x-axis limits
x_lim <- range(
     c(yb, theta_hat),
     na.rm = TRUE
)

x_pad <- 0.05 * diff(x_lim)

plot(
     x        = NA,
     y        = NA,
     xlim     = c(x_lim[1] - x_pad, x_lim[2] + x_pad),
     ylim     = c(1, 4),
     xlab     = "Score",
     ylab     = "",
     main     = "",
     yaxt     = "n",
     cex.axis = 0.8
)

axis(
     side   = 2,
     at     = c(2, 3),
     labels = c(expression(hat(theta)[j]), expression(bar(y)[j])),
     las    = 1
)

abline(
     h   = c(2, 3),
     col = c(4, 2),
     lwd = 2
)

for (j in seq_len(m)) {
     segments(
          x0  = theta_hat[j],
          y0  = 2,
          x1  = yb[j],
          y1  = 3,
          col = adjustcolor("black", alpha.f = 0.45),
          lwd = 1
     )
}

dev.off()

# End --------------------------------------------------------------------------