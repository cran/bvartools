## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data, fig.align='center', fig.height=5, fig.width=4.5--------------------
library(bvartools)

data("e1")
e1 <- diff(log(e1))

plot(e1) # Plot the series

## -----------------------------------------------------------------------------
data <- gen_var(e1, p = 2, deterministic = "const")

y <- data$Y[, 1:73]
x <- data$Z[, 1:73]

## -----------------------------------------------------------------------------
A_freq <- tcrossprod(y, x) %*% solve(tcrossprod(x)) # Calculate estimates
round(A_freq, 3) # Round estimates and print

## -----------------------------------------------------------------------------
u_freq <- y - A_freq %*% x
u_sigma_freq <- tcrossprod(u_freq) / (ncol(y) - nrow(x))
round(u_sigma_freq * 10^4, 2)

## ----flat prior---------------------------------------------------------------
# Reset random number generator for reproducibility
set.seed(1234567)

iter <- 10000 # Number of iterations of the Gibbs sampler
burnin <- 5000 # Number of burn-in draws
store <- iter - burnin

t <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
m <- k * nrow(x) # Number of estimated coefficients

# Set (uninformative) priors
a_mu_prior <- matrix(0, m) # Vector of prior parameter means
a_v_i_prior <- diag(0, m) # Inverse of the prior covariance matrix

u_sigma_df_prior <- 0 # Prior degrees of freedom
u_sigma_scale_prior <- diag(0, k) # Prior covariance matrix
u_sigma_df_post <- t + u_sigma_df_prior # Posterior degrees of freedom

# Initial values
u_sigma <- diag(.00001, k)
u_sigma_i <- solve(u_sigma)

# Data containers for posterior draws
draws_a <- matrix(NA, m, store)
draws_sigma <- matrix(NA, k^2, store)

# Start Gibbs sampler
for (draw in 1:iter) {
  # Draw conditional mean parameters
  a <- post_normal(y, x, u_sigma_i, a_mu_prior, a_v_i_prior)
  
  # Draw variance-covariance matrix
  u <- y - matrix(a, k) %*% x # Obtain residuals
  u_sigma_scale_post <- solve(u_sigma_scale_prior + tcrossprod(u))
  u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
  u_sigma <- solve(u_sigma_i) # Invert Sigma_i to obtain Sigma
  
  # Store draws
  if (draw > burnin) {
    draws_a[, draw - burnin] <- a
    draws_sigma[, draw - burnin] <- u_sigma
  }
}

## -----------------------------------------------------------------------------
A <- rowMeans(draws_a) # Obtain means for every row
A <- matrix(A, k) # Transform mean vector into a matrix
A <- round(A, 3) # Round values
dimnames(A) <- list(dimnames(y)[[1]], dimnames(x)[[1]]) # Rename matrix dimensions

A # Print

## -----------------------------------------------------------------------------
Sigma <- rowMeans(draws_sigma) # Obtain means for every row
Sigma <- matrix(Sigma, k) # Transform mean vector into a matrix
Sigma <- round(Sigma * 10^4, 2) # Round values
dimnames(Sigma) <- list(dimnames(y)[[1]], dimnames(y)[[1]]) # Rename matrix dimensions

Sigma # Print

## ----bvar-object--------------------------------------------------------------
bvar_est <- bvar(y = y, x = x, A = draws_a[1:18,],
                 C = draws_a[19:21, ], Sigma = draws_sigma)

## ----thin---------------------------------------------------------------------
bvar_est <- thin(bvar_est, thin = 5)

## ----forecasts, fig.width=5.5, fig.height=5.5---------------------------------
bvar_pred <- predict(bvar_est, n.ahead = 10, new_D = rep(1, 10))

plot(bvar_pred)

## ----feir, fig.width=5.5, fig.height=4.5--------------------------------------
FEIR <- irf(bvar_est, impulse = "income", response = "cons", n.ahead = 8)

plot(FEIR, main = "Forecast Error Impulse Response", xlab = "Period", ylab = "Response")

## ----oir, fig.width=5.5, fig.height=4.5---------------------------------------
OIR <- irf(bvar_est, impulse = "income", response = "cons", n.ahead = 8, type = "oir")

plot(OIR, main = "Orthogonalised Impulse Response", xlab = "Period", ylab = "Response")

## ----gir, fig.width=5.5, fig.height=4.5---------------------------------------
GIR <- irf(bvar_est, impulse = "income", response = "cons", n.ahead = 8, type = "gir")

plot(GIR, main = "Generalised Impulse Response", xlab = "Period", ylab = "Response")

## ----fevd-oir, fig.width=5.5, fig.height=4.5----------------------------------
bvar_fevd_oir <- fevd(bvar_est, response = "cons")

plot(bvar_fevd_oir, main = "OIR-based FEVD of consumption")

## ----fevd-gir, fig.width=5.5, fig.height=4.5----------------------------------
bvar_fevd_gir <- fevd(bvar_est, response = "cons", type = "gir")

plot(bvar_fevd_gir, main = "GIR-based FEVD of consumption")

