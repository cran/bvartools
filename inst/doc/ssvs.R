## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data, fig.align='center', fig.height=5, fig.width=4.5---------------
# devtools::install_github("franzmohr/bvartools")
library(bvartools)

# Load and transform data
data("e1")
e1 <- diff(log(e1))

# Generate VAR
data <- gen_var(e1, p = 4, deterministic = "const")

# Get data matrices
y <- data$Y[, 1:71]
x <- data$Z[, 1:71]

## ------------------------------------------------------------------------
# Reset random number generator for reproducibility
set.seed(1234567)

t <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
m <- k * nrow(x) # Number of estimated coefficients

# Coefficient priors
a_mu_prior <- matrix(0, m) # Vector of prior means

# SSVS priors (semiautomatic approach)
vs_prior <- ssvs_prior(data, semiautomatic = c(.1, 10))
tau0 <- vs_prior$tau0
tau1 <- vs_prior$tau1

# Prior for inclusion parameter
prob_prior <- matrix(0.5, m)

# Prior for variance-covariance matrix
u_sigma_df_prior <- 0 # Prior degrees of freedom
u_sigma_scale_prior <- diag(0, k) # Prior covariance matrix
u_sigma_df_post <- t + u_sigma_df_prior # Posterior degrees of freedom

## ------------------------------------------------------------------------
# Initial values
a <- matrix(0, m)
a_v_i_prior <- diag(1 / c(tau1)^2, m) # Inverse of the prior covariance matrix

# Data containers for posterior draws
iter <- 15000 # Number of total Gibs sampler draws
burnin <- 5000 # Number of burn-in draws

store <- iter - burnin
draws_a <- matrix(NA, m, store)
draws_lambda <- matrix(NA, m, store)
draws_sigma <- matrix(NA, k^2, store)

## ------------------------------------------------------------------------
# Reset random number generator for reproducibility
set.seed(1234567)

# Start Gibbs sampler
for (draw in 1:iter) {
  # Draw variance-covariance matrix
  u <- y - matrix(a, k) %*% x # Obtain residuals
  # Scale posterior
  u_sigma_scale_post <- solve(u_sigma_scale_prior + tcrossprod(u))
  # Draw posterior of inverse sigma
  u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
  # Obtain sigma
  u_sigma <- solve(u_sigma_i)
  
  # Draw conditional mean parameters
  a <- post_normal(y, x, u_sigma_i, a_mu_prior, a_v_i_prior)
  
  # Draw inclusion parameters and update priors
  temp <- ssvs(a, tau0, tau1, prob_prior, include = 1:36)
  a_v_i_prior <- temp$V_i # Update prior
  
  # Store draws
  if (draw > burnin) {
    draws_a[, draw - burnin] <- a
    draws_lambda[, draw - burnin] <- temp$lambda
    draws_sigma[, draw - burnin] <- u_sigma
  }
}

## ------------------------------------------------------------------------
A <- rowMeans(draws_a) # Obtain means for every parameter
A <- matrix(A, k) # Transform mean vector into matrix
A <- round(A, 3) # Round values
dimnames(A) <- list(dimnames(y)[[1]], dimnames(x)[[1]]) # Rename matrix dimensions

t(A) # Print

## ------------------------------------------------------------------------
lambda <- rowMeans(draws_lambda) # Obtain means for every row
lambda <- matrix(lambda, k) # Transform mean vector into a matrix
lambda <- round(lambda, 2) # Round values
dimnames(lambda) <- list(dimnames(y)[[1]], dimnames(x)[[1]]) # Rename matrix dimensions

t(lambda) # Print

## ---- fig.height=3.5, fig.width=4.5--------------------------------------
hist(draws_a[6,], main = "Consumption ~ First lag of income", xlab = "Value of posterior draw")

## ------------------------------------------------------------------------
# Select variables that should be included
include_var <- c(lambda >= .4)

# Update prior variances
diag(a_v_i_prior)[!include_var] <- 100000 # Very tight prior close to zero
diag(a_v_i_prior)[include_var] <- 1 / 9 # Relatively uninformative prior

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

## ------------------------------------------------------------------------
A <- rowMeans(draws_a) # Obtain means for every row
A <- matrix(A, k) # Transform mean vector into a matrix
A <- round(A, 3) # Round values
dimnames(A) <- list(dimnames(y)[[1]], dimnames(x)[[1]]) # Rename matrix dimensions

t(A) # Print

## ----bvar-object---------------------------------------------------------
bvar_est <- bvar(y = y, x = x, A = draws_a[1:36,],
                 C = draws_a[37:39, ], Sigma = draws_sigma)

## ----thin----------------------------------------------------------------
bvar_est <- thin(bvar_est, thin = 5)

## ----forecasts, fig.width=5.5, fig.height=5.5----------------------------
bvar_pred <- predict(bvar_est, n.ahead = 10, new_D = rep(1, 10))

plot(bvar_pred)

## ----oir, fig.width=5.5, fig.height=4.5----------------------------------
OIR <- irf(bvar_est, impulse = "income", response = "cons", n.ahead = 8, type = "oir")

plot(OIR, main = "Orthogonalised Impulse Response", xlab = "Period", ylab = "Response")

