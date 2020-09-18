## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data, fig.align='center', fig.height=5, fig.width=4.5--------------------
library(bvartools)

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Reduce number of oberservations
e1 <- window(e1, end = c(1978, 4))

# Plot the series
plot(e1)

## -----------------------------------------------------------------------------
model <- gen_var(e1, p = 2, deterministic = "const",
                 iterations = 5000, burnin = 1000)

## -----------------------------------------------------------------------------
model_with_priors <- add_priors(model,
                                coef = list(v_i = 0, v_i_det = 0),
                                sigma = list(df = 1, scale = .0001))

## ---- message=FALSE, warning=FALSE--------------------------------------------
bvar_est <- draw_posterior(model_with_priors)

## -----------------------------------------------------------------------------
summary(bvar_est)

## -----------------------------------------------------------------------------
# Obtain data for LS estimator
y <- t(model$data$Y)
z <- t(model$data$Z)

# Calculate LS estimates
A_freq <- tcrossprod(y, z) %*% solve(tcrossprod(z))

# Round estimates and print
round(A_freq, 3)

## -----------------------------------------------------------------------------
bvar_est <- thin_posterior(bvar_est, thin = 10)

## ----forecasts, fig.width=5.5, fig.height=5.5---------------------------------
bvar_pred <- predict(bvar_est, n.ahead = 10, new_d = rep(1, 10))

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

## ----flat prior---------------------------------------------------------------
# Reset random number generator for reproducibility
set.seed(1234567)

# Get data matrices
y <- t(model_with_priors$data$Y)
x <- t(model_with_priors$data$Z)

tt <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
m <- k * nrow(x) # Number of estimated coefficients

# Priors for coefficients
a_mu_prior <- model_with_priors$priors$coefficients$mu # Prior means
a_v_i_prior <- model_with_priors$priors$coefficients$v_i # Prior precisions

# Priors for error variance-covariance matrix
u_sigma_df_prior <- model_with_priors$priors$sigma$df # Prior degrees of freedom
u_sigma_scale_prior <- model_with_priors$priors$sigma$scale # Prior covariance matrix
u_sigma_df_post <- tt + u_sigma_df_prior # Posterior degrees of freedom

# Initial values for variance-covariance matrix
u_sigma <- diag(.00001, k)
u_sigma_i <- solve(u_sigma)

# Number of iterations of the Gibbs sampler
iterations <- model_with_priors$model$iterations 
# Number of burn-in draws
burnin <- model_with_priors$model$burnin
# Total number of draws
draws <- iterations + burnin

# Storate for posterior draws
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
  
  # Store draws
  if (draw > burnin) {
    draws_a[, draw - burnin] <- a
    draws_sigma[, draw - burnin] <- solve(u_sigma_i) # Invert Sigma_i to obtain Sigma
  }
}

## ----bvar-object--------------------------------------------------------------
bvar_est_two <- bvar(y = model_with_priors$data$Y,
                     x = model_with_priors$data$Z,
                     A = draws_a[1:18,],
                     C = draws_a[19:21, ],
                     Sigma = draws_sigma)

## -----------------------------------------------------------------------------
summary(bvar_est_two)

