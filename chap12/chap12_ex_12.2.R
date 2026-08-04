# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
setwd("~/Dropbox/UN/bayes_book_en")

# Load the required package
suppressMessages(suppressWarnings(library(gtools)))

# Parameters
alpha <- 10

H_values <- c(5, 10, 50)

G0 <- function(y) punif(y, min = 0, max = 1)

set.seed(123)

# Simulate a finite Dirichlet approximation for each value of H
for (H in H_values) {
     # Simulate the atoms
     vartheta <- runif(
          n   = H,
          min = 0,
          max = 1
     )
     
     # Simulate the weights
     omega <- as.numeric(
          gtools::rdirichlet(
               n     = 1,
               alpha = rep(alpha / H, H)
          )
     )
     
     # Sort the atoms and corresponding weights
     ord          <- order(vartheta)
     vartheta_ord <- vartheta[ord]
     omega_ord    <- omega[ord]
     
     # Construct the distribution function of G_H
     G_H <- stepfun(
          x     = vartheta_ord,
          y     = c(0, cumsum(omega_ord)),
          right = FALSE
     )
     
     # Create a separate file for each value of H
     pdf(
          file = paste0(
               "dp_aproximacion_dirichlet_finita_H_",
               H,
               ".pdf"
          ),
          width     = 5,
          height    = 5,
          pointsize = 18
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     # Plot the realization of G_H
     plot(
          G_H,
          xlim      = c(0, 1),
          ylim      = c(0, 1),
          do.points = FALSE,
          verticals = TRUE,
          xlab      = "y",
          ylab      = expression(G[H](y)),
          main      = "",
          lwd       = 1,
          col       = 4
     )
     
     # Add the base distribution function
     curve(
          G0,
          from = 0,
          to   = 1,
          n    = 1000,
          lwd  = 2,
          add  = TRUE
     )
     
     dev.off()
}

# End --------------------------------------------------------------------------