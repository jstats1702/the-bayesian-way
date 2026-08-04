# Settings ---------------------------------------------------------------------

rm(list = ls())

# Working directory
setwd("~/Dropbox/UN/bayes_book_en")

# Load the required package
suppressMessages(suppressWarnings(library(gtools)))

# Parameters
k <- 10000  # Number of points

alpha_values <- c(0.1, 1, 10, 100)  # Different values of alpha

G0 <- function(x) pnorm(x)  # Base measure: standard Normal distribution function

set.seed(123)

# Simulation for different values of alpha
for (alpha in alpha_values) {
     pdf(
          file      = paste0("dp_simulacion_particion_alpha_", alpha, ".pdf"),
          width     = 5,
          height    = 5,
          pointsize = 15
     )
     
     par(
          mfrow = c(1, 1),
          mar   = c(3, 3, 1.4, 1.4),
          mgp   = c(1.75, 0.75, 0)
     )
     
     # Configure the plot
     plot(
          NA,
          NA,
          xlim = c(-3, 3),
          ylim = c(0, 1),
          xlab = "x",
          ylab = "G(x)",
          main = ""
     )
     
     # Generate multiple realizations
     for (l in 1:10) {
          # Simulate and sort the x values
          x <- sort(
               runif(
                    n   = k,
                    min = -3,
                    max = 3
               )
          )
          
          # Compute the concentration parameters
          a <- numeric(k + 1)
          
          a[1]     <- alpha * G0(x[1])
          a[k + 1] <- alpha * (1 - G0(x[k]))
          
          for (j in 2:k) {
               a[j] <- alpha * (G0(x[j]) - G0(x[j - 1]))
          }
          
          # Simulate from the Dirichlet distribution
          u <- gtools::rdirichlet(
               n     = 1,
               alpha = a
          )
          
          # Plot the cumulative sum of the simulated weights
          lines(
               x,
               cumsum(u)[-length(u)],
               type = "l",
               col  = which(alpha_values == alpha)
          )
     }
     
     # Add the base measure curve
     curve(
          G0,
          from = -3,
          to   = 3,
          n    = 1000,
          lwd  = 2,
          add  = TRUE
     )
     
     dev.off()
}

# End --------------------------------------------------------------------------