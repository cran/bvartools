## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data, fig.align='center', fig.height=5, fig.width=4.5--------------------
library(bvartools)

# Load and transform data
data("e1")
e1 <- diff(log(e1))

# Shorten time series
e1 <- window(e1, end = c(1978, 4))

# Generate VAR
data <- gen_var(e1, p = 4, deterministic = "const",
                iterations = 10000, burnin = 5000)

## -----------------------------------------------------------------------------
# Reset random number generator for reproducibility
set.seed(1234567)

# Get data matrices
y <- t(data$data$Y)
x <- t(data$data$Z)

tt <- ncol(y) # Number of observations
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
u_sigma_scale_prior <- diag(0.00001, k) # Prior covariance matrix
u_sigma_df_post <- tt + u_sigma_df_prior # Posterior degrees of freedom

## -----------------------------------------------------------------------------
# Initial values
a <- matrix(0, m)
a_v_i_prior <- diag(1 / c(tau1)^2, m) # Inverse of the prior covariance matrix

# Data containers for posterior draws
iterations <- 10000 # Number of total Gibs sampler draws
burnin <- 5000 # Number of burn-in draws
draws <- iterations + burnin # Total number of draws

draws_a <- matrix(NA, m, iterations)
draws_lambda <- matrix(NA, m, iterations)
draws_sigma <- matrix(NA, k^2, iterations)

## -----------------------------------------------------------------------------
# Start Gibbs sampler
for (draw in 1:draws) {
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
  a_v_i_prior <- temp$v_i # Update prior
  
  # Store draws
  if (draw > burnin) {
    draws_a[, draw - burnin] <- a
    draws_lambda[, draw - burnin] <- temp$lambda
    draws_sigma[, draw - burnin] <- u_sigma
  }
}

## -----------------------------------------------------------------------------
bvar_est <- bvar(y = data$data$Y, x = data$data$Z,
                 A = list(coeffs = draws_a[1:36,],
                          lambda = draws_lambda[1:36,]),
                 C = list(coeffs = draws_a[37:39, ],
                          lambda = draws_lambda[37:39,]),
                 Sigma = draws_sigma)

bvar_summary <- summary(bvar_est)

bvar_summary

## ---- fig.height=3.5, fig.width=4.5-------------------------------------------
hist(draws_a[6,], main = "Consumption ~ First lag of income", xlab = "Value of posterior draw")

## -----------------------------------------------------------------------------
# Get inclusion probabilities
lambda <- bvar_summary$coefficients$lambda

# Select variables that should be included
include_var <- c(lambda >= .4)

# Update prior variances
diag(a_v_i_prior)[!include_var] <- 1 / 0.00001 # Very tight prior close to zero
diag(a_v_i_prior)[include_var] <- 1 / 9 # Relatively uninformative prior

# Data containers for posterior draws
draws_a <- matrix(NA, m, iterations)
draws_sigma <- matrix(NA, k^2, iterations)

# Start Gibbs sampler
for (draw in 1:draws) {
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
bvar_est <- bvar(y = data$data$Y, x = data$data$Z, A = draws_a[1:36,],
                 C = draws_a[37:39, ], Sigma = draws_sigma)

summary(bvar_est)

## ---- eval = FALSE------------------------------------------------------------
#  # Obtain priors
#  model_with_priors <- add_priors(data,
#                                  ssvs = list(inprior = 0.5, semiautomatic = c(0.01, 10), exclude_det = TRUE),
#                                  sigma = list(df = 0, scale = 0.00001))

## ---- message = FALSE, warning = FALSE, eval = FALSE--------------------------
#  ssvs_est <- draw_posterior(model_with_priors)

