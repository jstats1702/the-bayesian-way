# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
setwd("~/Dropbox/UN/bayes_book_en")

# Packages
library(maps)
library(sf)
library(spdep)
library(Matrix)
library(ggplot2)


# Data -------------------------------------------------------------------------

flu <- read.csv("flu_data.csv", header = TRUE)    # Complete dataset
flu <- flu[485:496, ]                             # First 12 weeks of 2013
flu <- flu[, -c(1, 2, 12)]                        # Excluding date, Alaska, and Hawaii

states.names <- colnames(flu)

data <- as.matrix(flu)
p    <- nrow(data)                                # T: number of weeks
m    <- ncol(data)                                # I: number of states


# Regions of the United States Department of Health and Human Services ---------

regions <- read.csv("flu_regions.csv", header = TRUE)     # Complete dataset
regions <- regions[-c(2, 12), ]                           # Excluding Alaska and Hawaii
regions <- as.matrix(regions)

e <- length(unique(regions))                              # Number of regions


# Distribution by race ---------------------------------------------------------

whites <- read.csv("flu_race.csv", header = TRUE)    # Complete dataset
whites <- whites[, 2]                                # Percentage of White population
whites <- whites[-c(2, 12)]                          # Excluding Alaska and Hawaii
whites <- as.matrix(whites)


# Distribution by age ----------------------------------------------------------

old <- read.csv("flu_age.csv", header = TRUE)    # Complete dataset
old <- old[, 4]                                  # Percentage of population older than 65 years
old <- old[-c(2, 12)]                            # Excluding Alaska and Hawaii
old <- as.matrix(old)


# Binary response --------------------------------------------------------------

y <- matrix(0, p, m)

for (i in 1:p) {
     for (j in 1:m) {
          if (data[i, j] > 7500) {
               y[i, j] <- 1
          }
     }
}

y <- as.matrix(c(y))


# Neighborhood structure among states -----------------------------------------

# Use planar geometry instead of s2
sf::sf_use_s2(FALSE)

# Load state polygons
usa.state <- map(
     database = "state",
     fill     = TRUE,
     plot     = FALSE
)

# Extract state identifiers
state.ID <- sapply(
     strsplit(usa.state$names, ":"),
     function(x) x[1]
)

# Convert the map to an sf object
usa.sf <- sf::st_as_sf(
     usa.state,
     IDs = state.ID
)

# Correct invalid geometries
usa.sf <- sf::st_make_valid(usa.sf)

# Project to planar coordinates
usa.sf <- sf::st_transform(
     usa.sf,
     crs = 5070
)

# Construct neighbor list
usa.nb <- spdep::poly2nb(
     usa.sf,
     queen = TRUE
)

# Convert neighborhoods to a binary matrix
usa.adj.mat <- A <- spdep::nb2mat(
     usa.nb,
     style = "B"
)   # Neighborhood matrix


# Neighborhood structure among regions ----------------------------------------

A2 <- c(
     0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
     1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     0, 1, 0, 1, 1, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 1, 1, 1, 0, 0, 0,
     0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
     0, 0, 0, 1, 0, 0, 1, 1, 1, 0,
     0, 0, 0, 1, 1, 1, 0, 1, 0, 0,
     0, 0, 0, 0, 1, 1, 1, 0, 1, 1,
     0, 0, 0, 0, 0, 1, 0, 1, 0, 1,
     0, 0, 0, 0, 0, 0, 0, 1, 1, 0
)

A2 <- matrix(
     A2,
     nrow  = e,
     ncol  = e,
     byrow = TRUE
)


# Plot of the neighborhood matrix among states ---------------------------------

pdf(
     file      = "flu_matriz_vecindad_estados.pdf",
     height    = 5,
     width     = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(1.4, 1.4, 0.4, 0.4),
     mgp   = c(1.75, 0.75, 0)
)

colorscale <- c("white", rev(heat.colors(100)))

image(
     x    = 1:m,
     y    = 1:m,
     z    = A,
     axes = FALSE,
     xlab = "",
     ylab = "",
     main = "",
     col  = colorscale[seq(floor(100 * min(A)), floor(100 * max(A)))]
)

axis(
     side     = 1,
     at       = 1:m,
     labels   = NA,
     las      = 2,
     cex.axis = 0.4
)

axis(
     side     = 2,
     at       = 1:m,
     labels   = NA,
     las      = 2,
     cex.axis = 0.4
)

box()

dev.off()


# Plot of the neighborhood matrix among regions --------------------------------

pdf(
     file      = "flu_matriz_vecindad_regiones.pdf",
     height    = 5,
     width     = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(1.4, 1.4, 0.4, 0.4),
     mgp   = c(1.75, 0.75, 0)
)

colorscale <- c("white", rev(heat.colors(100)))

image(
     x    = 1:e,
     y    = 1:e,
     z    = A2,
     axes = FALSE,
     xlab = "",
     ylab = "",
     main = "",
     col  = colorscale[seq(floor(100 * min(A2)), floor(100 * max(A2)))]
)

axis(
     side     = 1,
     at       = 1:e,
     labels   = 1:e,
     las      = 2,
     cex.axis = 0.75
)

axis(
     side     = 2,
     at       = 1:e,
     labels   = 1:e,
     las      = 2,
     cex.axis = 0.75
)

box()

dev.off()


# State precision matrix -------------------------------------------------------

W <- -A

diag(W) <- rowSums(A)

r     <- rankMatrix(W)[1]
n.nei <- rowSums(A)


# Region precision matrix ------------------------------------------------------

W2 <- -A2

diag(W2) <- rowSums(A2)

r2     <- rankMatrix(W2)[1]
n.nei2 <- rowSums(A2)


# Map of the proportion of White population -----------------------------------

# Base state map
all_states <- map_data("state")

# Merge covariate with map
whites <- read.csv("flu_race.csv", header = TRUE)    # Complete dataset
whites <- whites[, 2]                                # Percentage of White population
whites <- whites[-c(2, 12)]                          # Excluding Alaska and Hawaii
whites <- as.data.frame(whites)
whites <- cbind(rownames(A), whites)

colnames(whites) <- c("region", "percent")

whites <- merge(
     all_states,
     whites,
     by = "region"
)

head(whites)

# Visualization
q <- ggplot() +
     geom_polygon(
          data   = whites,
          aes(
               x     = long,
               y     = lat,
               group = group,
               fill  = percent
          ),
          colour = "black"
     ) +
     scale_fill_continuous(
          low   = "thistle2",
          high  = "darkred",
          guide = "colorbar"
     ) +
     labs(
          fill  = "%",
          title = "",
          x     = "",
          y     = ""
     ) +
     theme(
          panel.background = element_rect(colour = "black"),
          panel.grid.major = element_line(size = 1),
          plot.margin      = margin(0, 0, 0, 0)
     ) +
     theme_minimal(base_size = 30)

q

ggsave(
     filename = "flu_mapa_proporcion_blancos.pdf",
     plot     = q,
     width    = 10,
     height   = 7
)

dev.off()


# Map of the proportion of population older than 65 years ---------------------

# Base state map
all_states <- map_data("state")

# Merge covariate with map
old <- read.csv("flu_age.csv", header = TRUE)    # Complete dataset
old <- old[, 4]                                  # Percentage of population older than 65 years
old <- old[-c(2, 12)]                            # Excluding Alaska and Hawaii
old <- as.data.frame(old)
old <- cbind(rownames(A), old)

colnames(old) <- c("region", "percent")

old <- merge(
     all_states,
     old,
     by = "region"
)

head(old)

# Visualization
q <- ggplot() +
     geom_polygon(
          data   = old,
          aes(
               x     = long,
               y     = lat,
               group = group,
               fill  = percent
          ),
          colour = "black"
     ) +
     scale_fill_continuous(
          low   = "yellow1",
          high  = "yellow4",
          guide = "colorbar"
     ) +
     labs(
          fill  = "%",
          title = "",
          x     = "",
          y     = ""
     ) +
     theme(
          panel.background = element_rect(colour = "black"),
          panel.grid.major = element_line(size = 1),
          plot.margin      = margin(0, 0, 0, 0)
     ) +
     theme_minimal(base_size = 30)

q

ggsave(
     filename = "flu_mapa_proporcion_mayores_65.pdf",
     plot     = q,
     width    = 10,
     height   = 7
)

dev.off()


# Map of the response by week --------------------------------------------------

# Base state map
all_states <- map_data("state")

# Binary response
y <- matrix(0, p, m)

for (i in 1:p) {
     for (j in 1:m) {
          if (data[i, j] > 7500) {
               y[i, j] <- 1
          }
     }
}

y <- t(y)

# Visualization by week
for (tt in 1:p) {
     
     # Merge response with map
     yt <- as.data.frame(as.factor(y[, tt]))
     yt <- cbind(rownames(A), yt)
     
     colnames(yt) <- c("region", paste("week_", tt, sep = ""))
     
     yt <- merge(
          all_states,
          yt,
          by = "region"
     )
     
     q <- ggplot() +
          geom_polygon(
               data   = yt,
               aes(
                    x     = long,
                    y     = lat,
                    group = group,
                    fill  = yt[, 7]
               ),
               colour = "black"
          ) +
          scale_fill_manual(
               values = c("white", "red"),
               guide  = "none"
          ) +
          labs(
               title = "",
               x     = "",
               y     = ""
          ) +
          theme_minimal(base_size = 30) +
          theme(
               legend.position  = "none",
               axis.text        = element_blank(),
               axis.ticks       = element_blank(),
               panel.grid       = element_blank(),
               panel.background = element_rect(colour = "black"),
               plot.margin      = margin(0, 30, 0, 0)
          )
     
     q
     
     ggsave(
          filename = paste("flu_mapa_respuesta_semana_", tt, ".pdf", sep = ""),
          plot     = q,
          width    = 10,
          height   = 7
     )
}

# End --------------------------------------------------------------------------

# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
setwd("~/Dropbox/UN/bayes_book_en")

# Packages
library(maps)
library(sf)
library(spdep)
library(Matrix)
library(ggplot2)


# Data -------------------------------------------------------------------------

flu <- read.csv("flu_data.csv", header = TRUE)    # Complete dataset
flu <- flu[485:496, ]                             # First 12 weeks of 2013
flu <- flu[, -c(1, 2, 12)]                        # Excluding date, Alaska, and Hawaii

states.names <- colnames(flu)

data <- as.matrix(flu)
p    <- nrow(data)                                # T: number of weeks
m    <- ncol(data)                                # I: number of states


# Regions of the United States Department of Health and Human Services ---------

regions <- read.csv("flu_regions.csv", header = TRUE)     # Complete dataset
regions <- regions[-c(2, 12), ]                           # Excluding Alaska and Hawaii
regions <- as.matrix(regions)

e <- length(unique(regions))                               # Number of regions


# Distribution by race ---------------------------------------------------------

whites <- read.csv("flu_race.csv", header = TRUE)    # Complete dataset
whites <- whites[, 2]                                # Percentage of White population
whites <- whites[-c(2, 12)]                          # Excluding Alaska and Hawaii
whites <- as.matrix(whites)


# Distribution by age ----------------------------------------------------------

old <- read.csv("flu_age.csv", header = TRUE)    # Complete dataset
old <- old[, 4]                                  # Percentage of population older than 65 years
old <- old[-c(2, 12)]                            # Excluding Alaska and Hawaii
old <- as.matrix(old)


# Binary response --------------------------------------------------------------

y <- matrix(0, p, m)

for (i in 1:p) {
     for (j in 1:m) {
          if (data[i, j] > 7500) {
               y[i, j] <- 1
          }
     }
}

y <- as.matrix(c(y))


# State neighborhood matrix ----------------------------------------------------

# Base state map
usa.state <- maps::map(
     database = "state",
     fill     = TRUE,
     plot     = FALSE
)

# Convert to an sf object
usa.sf <- sf::st_as_sf(usa.state)

# State identifier
usa.sf$region <- usa.sf$ID

# Correct geometries
usa.sf <- sf::st_make_valid(usa.sf)

# Temporarily disable the s2 engine
s2.original <- sf::sf_use_s2()
sf::sf_use_s2(FALSE)

# Neighborhood structure among states
usa.nb <- spdep::poly2nb(
     usa.sf,
     row.names = usa.sf$region,
     queen     = TRUE
)

# Adjacency matrix
usa.adj.mat <- A <- spdep::nb2mat(
     usa.nb,
     style       = "B",
     zero.policy = TRUE
)

rownames(A) <- usa.sf$region
colnames(A) <- usa.sf$region

# Restore the original configuration
sf::sf_use_s2(s2.original)


# Region neighborhood matrix --------------------------------------------------

A2 <- matrix(
     c(
          0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
          1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
          0, 1, 0, 1, 1, 0, 0, 0, 0, 0,
          0, 0, 1, 0, 1, 1, 1, 0, 0, 0,
          0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
          0, 0, 0, 1, 0, 0, 1, 1, 1, 0,
          0, 0, 0, 1, 1, 1, 0, 1, 0, 0,
          0, 0, 0, 0, 1, 1, 1, 0, 1, 1,
          0, 0, 0, 0, 0, 1, 0, 1, 0, 1,
          0, 0, 0, 0, 0, 0, 0, 1, 1, 0
     ),
     nrow  = e,
     ncol  = e,
     byrow = TRUE
)


# Precision matrix W for states -----------------------------------------------

W <- -A
diag(W) <- rowSums(A)

r     <- rankMatrix(W)[1]
n.nei <- rowSums(A)


# Precision matrix W2 for regions ---------------------------------------------

W2 <- -A2
diag(W2) <- rowSums(A2)

r2     <- rankMatrix(W2)[1]
n.nei2 <- rowSums(A2)


# Design matrix X for fixed effects --------------------------------------------

k <- 4

X <- matrix(0, m * p, k)

X[, 1] <- 1                                           # Intercept
X[, 2] <- whites %x% as.matrix(rep(1, p))             # Percentage of White population
X[, 3] <- old %x% as.matrix(rep(1, p))                # Percentage of older population
X[, 4] <- as.matrix(rep(1, m)) %x% as.matrix(1:p)     # Week

M <- dim(X)[1]
k <- dim(X)[2]


# Design matrix Z for spatial effects ------------------------------------------

k2 <- e

Z <- matrix(0, m * p, k2)

for (i in 1:m) {
     Z[((i - 1) * p + 1):(i * p), regions[i]] <- 1
}

rm(i)

M  <- dim(Z)[1]
k2 <- dim(Z)[2]

head(Z, 2 * p)


# Gibbs sampler ----------------------------------------------------------------

# Implementation note
# The data are stored by state, allowing weeks to vary within each state;
# for this reason, the temporal precisions are organized as rep(ka.w, m).
# In addition, although xi is theoretically presented as a separate parameter,
# it is incorporated in the code as the fourth column of X and updated jointly with beta.

sample.w <- function(X, Z, beta, gama, theta, phi, ka.w, y)
{
     a <- b <- rep(0, M)
     
     a[y == 0] <- -Inf
     b[y == 1] <- Inf
     
     mu <- X %*% beta + Z %*% gama + theta + phi
     w  <- matrix(NA, M, 1)
     
     for (tt in 1:p) {
          idx <- tt + (0:(m - 1)) * p
          
          w[idx] <- truncnorm::rtruncnorm(
               length(idx),
               a[idx],
               b[idx],
               mu[idx],
               sqrt(1 / ka.w[tt])
          )
     }
     
     as.numeric(w)
}


sample.beta <- function(X, Z, sig2.0, w, gama, theta, phi, ka.w)
{
     V.inv <- diag(rep(ka.w, m))
     
     SIGMA <- chol2inv(
          chol(
               diag(k) / sig2.0 +
                    t(X) %*% V.inv %*% X
          )
     )
     
     MEAN <- SIGMA %*% t(X) %*% V.inv %*% (
          w - Z %*% gama - theta - phi
     )
     
     as.numeric(mvtnorm::rmvnorm(1, MEAN, SIGMA))
}


sample.gama <- function(X, Z, W2, sig2.0, w, beta, theta, phi, ka.w)
{
     V.inv <- diag(rep(ka.w, m))
     
     SIGMA <- chol2inv(
          chol(
               W2 / sig2.0 +
                    t(Z) %*% V.inv %*% Z
          )
     )
     
     MEAN <- SIGMA %*% t(Z) %*% V.inv %*% (
          w - X %*% beta - theta - phi
     )
     
     gama <- as.numeric(mvtnorm::rmvnorm(1, MEAN, SIGMA))
     gama <- gama - mean(gama)
     
     gama
}


sample.theta <- function(X, Z, beta, gama, w, phi, ka.h, ka.w)
{
     theta <- matrix(NA, M, 1)
     
     for (tt in 1:p) {
          idx <- tt + (0:(m - 1)) * p
          
          SIGMA <- diag(m) / (ka.h[tt] + ka.w[tt])
          
          MEAN <- SIGMA %*% (
               ka.w[tt] * (
                    w[idx] -
                         X[idx, ] %*% beta -
                         Z[idx, ] %*% gama -
                         phi[idx]
               )
          )
          
          theta[idx] <- as.numeric(mvtnorm::rmvnorm(1, MEAN, SIGMA))
          # theta[idx] <- theta[idx] - mean(theta[idx])
     }
     
     as.numeric(theta)
}


sample.phi <- function(X, Z, W, beta, gama, w, theta, ka.c, ka.w)
{
     phi <- matrix(NA, M, 1)
     
     for (tt in 1:p) {
          idx <- tt + (0:(m - 1)) * p
          
          SIGMA <- chol2inv(
               chol(
                    ka.c[tt] * W +
                         ka.w[tt] * diag(m)
               )
          )
          
          MEAN <- SIGMA %*% (
               ka.w[tt] * (
                    w[idx] -
                         X[idx, ] %*% beta -
                         Z[idx, ] %*% gama -
                         theta[idx]
               )
          )
          
          phi[idx] <- as.numeric(mvtnorm::rmvnorm(1, MEAN, SIGMA))
          phi[idx] <- phi[idx] - mean(phi[idx])
     }
     
     as.numeric(phi)
}


sample.ka.w <- function(nu.0, X, Z, beta, gama, w, theta, phi)
{
     ka.w  <- matrix(NA, p, 1)
     shape <- nu.0 / 2 + m / 2
     
     for (tt in 1:p) {
          idx <- tt + (0:(m - 1)) * p
          
          rate <- nu.0 / 2 +
               0.5 * sum(
                    (
                         w[idx] -
                              X[idx, ] %*% beta -
                              Z[idx, ] %*% gama -
                              theta[idx] -
                              phi[idx]
                    )^2
               )
          
          ka.w[tt] <- rgamma(1, shape, rate)
     }
     
     as.numeric(ka.w)
}


sample.ka.h <- function(a.h, b.h, theta)
{
     ka.h  <- matrix(NA, p, 1)
     shape <- a.h + m / 2
     
     for (tt in 1:p) {
          idx <- tt + (0:(m - 1)) * p
          
          rate <- b.h + 0.5 * sum(theta[idx]^2)
          
          ka.h[tt] <- rgamma(1, shape, rate)
     }
     
     as.numeric(ka.h)
}


sample.ka.c <- function(a.c, b.c, W, phi)
{
     ka.c  <- matrix(NA, p, 1)
     shape <- a.c + r / 2
     
     for (tt in 1:p) {
          idx <- tt + (0:(m - 1)) * p
          
          rate <- b.c + 0.5 * as.numeric(
               t(phi[idx]) %*% W %*% phi[idx]
          )
          
          ka.c[tt] <- rgamma(1, shape, rate)
     }
     
     as.numeric(ka.c)
}

MCMC <- function(y, X, Z, W, W2, sig2.0, nu.0, a.h, b.h, a.c, b.c,
                 n.sams, n.burn, n.skip)
{
     require(Matrix)
     library(mvtnorm)
     
     # Dimensions
     m  <- dim(W)[1]          # Number of states
     M  <- dim(X)[1]          # Number of observations
     k  <- dim(X)[2]          # Number of fixed covariates
     k2 <- dim(Z)[2]          # Number of spatial covariates
     p  <- M / m              # Number of time points
     r  <- rankMatrix(W)[1]   # Rank of W, states
     r2 <- rankMatrix(W2)[1]  # Rank of W2, regions
     
     # Total number of iterations
     n.total <- n.burn + n.skip * n.sams
     
     # Storage of posterior samples
     beta.chain   <- matrix(NA, n.sams, k)
     gama.chain   <- matrix(NA, n.sams, k2)
     w.chain      <- matrix(NA, n.sams, M)
     theta.chain  <- matrix(NA, n.sams, M)
     phi.chain    <- matrix(NA, n.sams, M)
     ka.w.chain   <- matrix(NA, n.sams, p)
     ka.h.chain   <- matrix(NA, n.sams, p)
     ka.c.chain   <- matrix(NA, n.sams, p)
     loglik.chain <- rep(NA, n.sams)
     
     # Parameter initialization
     beta  <- rep(0, k)
     gama  <- rep(0, k2)
     w     <- rep(0, M)
     theta <- rep(0, M)
     phi   <- rep(0, M)
     ka.w  <- rep(1, p)
     ka.h  <- rep(1, p)
     ka.c  <- rep(1, p)
     
     # Progress points
     progreso <- unique(round(seq(0.10, 1, by = 0.10) * n.total))
     
     # MCMC
     s <- 0
     
     for (l in 1:n.total) {
          
          # Sampling
          w     <- sample.w(X, Z, beta, gama, theta, phi, ka.w, y)
          beta  <- sample.beta(X, Z, sig2.0, w, gama, theta, phi, ka.w)
          gama  <- sample.gama(X, Z, W2, sig2.0, w, beta, theta, phi, ka.w)
          theta <- sample.theta(X, Z, beta, gama, w, phi, ka.h, ka.w)
          phi   <- sample.phi(X, Z, W, beta, gama, w, theta, ka.c, ka.w)
          ka.w  <- sample.ka.w(nu.0, X, Z, beta, gama, w, theta, phi)
          ka.h  <- sample.ka.h(a.h, b.h, theta)
          ka.c  <- sample.ka.c(a.c, b.c, W, phi)
          
          # Posterior storage
          if (l > n.burn && (l - n.burn) %% n.skip == 0) {
               s <- s + 1
               
               beta.chain[s, ]  <- beta
               gama.chain[s, ]  <- gama
               w.chain[s, ]     <- w
               theta.chain[s, ] <- theta
               phi.chain[s, ]   <- phi
               ka.w.chain[s, ]  <- ka.w
               ka.h.chain[s, ]  <- ka.h
               ka.c.chain[s, ]  <- ka.c
               
               eta <- X %*% beta + Z %*% gama + theta + phi
               pi  <- pnorm(sqrt(rep(ka.w, m)) * eta)
               
               loglik.chain[s] <- sum(dbinom(y, size = 1, prob = pi, log = TRUE))
          }
          
          # Progress
          if (l %in% progreso || l == n.total) {
               cat(
                    "Progreso: ",
                    round(100 * l / n.total, 1),
                    "%\n",
                    sep = ""
               )
          }
     }
     
     # Output
     list(
          beta.chain   = beta.chain,
          gama.chain   = gama.chain,
          w.chain      = w.chain,
          theta.chain  = theta.chain,
          phi.chain    = phi.chain,
          ka.w.chain   = ka.w.chain,
          ka.h.chain   = ka.h.chain,
          ka.c.chain   = ka.c.chain,
          loglik.chain = loglik.chain,
          n.eff        = n.sams
     )
}


# Hyperparameters --------------------------------------------------------------

# Prior specification
sig2.0 <- 100  # psi0 = 1/sig2.0
nu.0   <- 3

# Target prior scale for the random effects on the latent scale
s0 <- 2

# Precision of theta_it
a.h <- 2
b.h <- (a.h - 1) * s0^2

# Precision of phi_it
a.c <- 2
b.c <- (0.7^2) * mean(rowSums(A)) * b.h  # BCG, p. 156


# Model fitting ----------------------------------------------------------------

n.sams <- 10000
n.burn <- 10000
n.skip <- 10

# Model fitting

set.seed(123)

fit <- MCMC(
     y, X, Z, W, W2,
     sig2.0, nu.0,
     a.h, b.h,
     a.c, b.c,
     n.sams, n.burn, n.skip
)

save(fit, file = "flu_muestras_distribucion_posterior.RData")

load(file = "flu_muestras_distribucion_posterior.RData")


# Log-likelihood trace ---------------------------------------------------------

plot(
     x    = fit$loglik.chain,
     type = "p",
     pch  = 16,
     cex  = 0.3,
     xlab = "Iteration",
     ylab = "Log-likelihood",
     main = "Log-likelihood"
)

dev.off()


# Posterior inference for the fixed effects and xi -----------------------------

resumen_posterior <- function(x)
{
     c(
          media = mean(x),
          sd    = sd(x),
          li95  = unname(quantile(x, 0.025)),
          ls95  = unname(quantile(x, 0.975))
     )
}

tabla.efectos <- t(
     apply(
          fit$beta.chain,
          2,
          resumen_posterior
     )
)

rownames(tabla.efectos) <- c(
     "beta_1",
     "beta_2",
     "beta_3",
     "xi"
)

tabla.efectos <- round(tabla.efectos, 3)

tabla.efectos


# Posterior inference for gamma -----------------------------------------------

# Posterior summary of the random effects
gama.media <- apply(fit$gama.chain, 2, mean)
gama.li95  <- apply(fit$gama.chain, 2, quantile, probs = 0.025)
gama.ls95  <- apply(fit$gama.chain, 2, quantile, probs = 0.975)

# Colors according to interval position
col.intervalo <- rep("black", length(gama.media))
col.intervalo[gama.li95 > 0] <- "darkgreen"
col.intervalo[gama.ls95 < 0] <- "red"

pdf(
     file      = "flu_inferencia_efectos_aleatorios.pdf",
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
     x    = 1:length(gama.media),
     y    = gama.media,
     type = "p",
     pch  = 16,
     col  = col.intervalo,
     ylim = range(c(gama.li95, gama.ls95)),
     xaxt = "n",
     xlab = "Random effect",
     ylab = "Posterior mean",
     main = ""
)

axis(
     side   = 1,
     at     = 1:length(gama.media),
     labels = 1:length(gama.media),
     las    = 1
)

arrows(
     x0     = 1:length(gama.media),
     y0     = gama.li95,
     x1     = 1:length(gama.media),
     y1     = gama.ls95,
     angle  = 90,
     code   = 3,
     length = 0,
     col    = col.intervalo
)

abline(
     h   = 0,
     lty = 2,
     lwd = 2
)

box()

dev.off()


# Posterior inference for the spatiotemporal effects ---------------------------

# Total spatiotemporal effect
efecto.et <- fit$theta.chain + fit$phi.chain

# Posterior mean by state and week
media.et <- matrix(NA, nrow = m, ncol = p)

for (i in 1:m) {
     idx <- ((i - 1) * p + 1):(i * p)
     media.et[i, ] <- apply(efecto.et[, idx], 2, mean)
}

# Region for each state
region.estado <- as.numeric(
     Z[seq(1, m * p, by = p), ] %*% (1:ncol(Z))
)

# Regions
regiones <- sort(unique(region.estado))

# Panels by region
for (rr in regiones) {
     
     idx.region <- which(region.estado == rr)
     
     pdf(
          file      = paste0("flu_inferencia_efectos_espaciotemporales_region_", rr, ".pdf"),
          width     = 7,
          height    = 5,
          pointsize = 25
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     plot(
          x    = 1:p,
          y    = media.et[idx.region[1], ],
          type = "n",
          ylim = range(media.et),
          xlab = "Week",
          ylab = "Effect",
          main = ""
     )
     
     for (i in idx.region) {
          lines(
               x   = 1:p,
               y   = media.et[i, ],
               col = 2,
               lwd = 1.5
          )
     }
     
     abline(
          h   = 0,
          lty = 2,
          lwd = 1.2
     )
     
     box()
     
     dev.off()
}


# Posterior inference for precision components ---------------------------------

# Posterior means of the precision components
ka.w.media <- apply(fit$ka.w.chain, 2, mean)
ka.h.media <- apply(fit$ka.h.chain, 2, mean)
ka.c.media <- apply(fit$ka.c.chain, 2, mean)

# Output file
pdf(
     file      = "flu_inferencia_componentes_precision.pdf",
     width     = 7,
     height    = 5,
     pointsize = 15
)

par(
     mfrow = c(1, 1),
     mar   = c(3, 3, 1.4, 1.4),
     mgp   = c(1.75, 0.75, 0)
)

# Plot
plot(
     x    = 1:p,
     y    = ka.w.media,
     type = "n",
     ylim = range(c(ka.w.media, ka.h.media, ka.c.media)),
     xaxt = "n",
     xlab = "Week",
     ylab = "Precision",
     main = ""
)

axis(
     side   = 1,
     at     = 1:p,
     labels = 1:p
)

lines(
     x   = 1:p,
     y   = ka.w.media,
     col = "black",
     lwd = 2
)

lines(
     x   = 1:p,
     y   = ka.h.media,
     col = "red",
     lwd = 2
)

lines(
     x   = 1:p,
     y   = ka.c.media,
     col = "blue",
     lwd = 2
)

legend(
     "topleft",
     legend = c(
          expression(kappa[t]),
          expression(tau[t]),
          expression(lambda[t])
     ),
     col = c("black", "red", "blue"),
     lty = 1,
     lwd = 2,
     bty = "n"
)

box()

dev.off()


# Prediction -------------------------------------------------------------------

post.pred <- function(i, tt, fit)
{
     n.sam    <- fit$n.eff
     beta.sam <- fit$beta.chain
     gama.sam <- fit$gama.chain
     w.it.sam <- rep(NA, n.sam)
     
     x.it    <- rep(0, k)
     x.it[1] <- 1
     x.it[2] <- whites[i]  # Percentage of White population
     x.it[3] <- old[i]     # Percentage of older population
     x.it[4] <- tt         # Week
     
     z.it <- rep(0, k2)
     z.it[regions[i]] <- 1
     
     idx <- (i - 1) * p + tt
     
     for (b in 1:n.sam) {
          kw.t     <- fit$ka.w.chain[b, tt]
          theta.it <- fit$theta.chain[b, idx]
          phi.it   <- fit$phi.chain[b, idx]
          
          mu <- sum(x.it * beta.sam[b, ]) +
               sum(z.it * gama.sam[b, ]) +
               theta.it +
               phi.it
          
          w.it.sam[b] <- rnorm(1, mu, sqrt(1 / kw.t))
     }
     
     y.it <- rep(0, n.sam)
     y.it[w.it.sam > 0] <- 1
     
     list(pr = mean(y.it), sam = y.it)
}

# Base state map
all_states <- map_data("state")

# Posterior predictive probability maps
for (tt in 1:p) {
     
     post.prob <- rep(NA, m)
     
     for (i in 1:m) {
          post.prob[i] <- post.pred(i, tt, fit)$pr
     }
     
     # Merge predictive probability with map
     post.prob <- as.data.frame(post.prob)
     post.prob <- cbind(rownames(A), post.prob)
     
     colnames(post.prob) <- c("region", "prob")
     
     post.prob <- merge(
          all_states,
          post.prob,
          by = "region"
     )
     
     q <- ggplot() +
          geom_polygon(
               data = post.prob,
               aes(
                    x     = long,
                    y     = lat,
                    group = group,
                    fill  = prob
               ),
               colour = "black"
          ) +
          scale_fill_continuous(
               low    = "white",
               high   = "red",
               limits = c(0, 1),
               guide  = "none"
          ) +
          labs(
               title = "",
               x     = "",
               y     = ""
          ) +
          theme_minimal(base_size = 30) +
          theme(
               legend.position  = "none",
               axis.text        = element_blank(),
               axis.ticks       = element_blank(),
               panel.grid       = element_blank(),
               panel.background = element_rect(colour = "black"),
               plot.margin      = margin(0, 30, 0, 0)
          )
     
     q
     
     ggsave(
          filename = paste("flu_inferencia_prediccion_semana_", tt, ".pdf", sep = ""),
          plot     = q,
          width    = 10,
          height   = 7
     )
     
     rm(post.prob)
}

# End --------------------------------------------------------------------------