## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data, fig.align='center', fig.height=5, fig.width=4.5--------------------
library(bvartools)

data("e6")

plot(e6) # Plot the series

## -----------------------------------------------------------------------------
data <- gen_vec(e6, p = 4, r = 1,
                const = "unrestricted", season = "unrestricted",
                iterations = 5000, burnin = 1000)

## -----------------------------------------------------------------------------
data <- add_priors(data,
                   coint = list(v_i = 0, p_tau_i = 1),
                   coef = list(v_i = 0, v_i_det = 0),
                   sigma = list(df = 0, scale = .0001))

## ----flat prior---------------------------------------------------------------
# Reset random number generator for reproducibility
set.seed(7654321)

# Obtain data matrices
y <- t(data$data$Y)
w <- t(data$data$W)
x <- t(data$data$X)

r <- data$model$rank # Set rank

tt <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
k_w <- nrow(w) # Number of regressors in error correction term
k_x <- nrow(x) # Number of differenced regressors and unrestrictec deterministic terms
k_gamma <- k * k_x # Total number of non-cointegration coefficients

k_alpha <- k * r # Number of elements in alpha
k_beta <- k_w * r # Number of elements in beta

# Priors
a_mu_prior <- data$priors$noncointegration$mu # Prior means
a_v_i_prior <- data$priors$noncointegration$v_i # Inverse of the prior covariance matrix

v_i <- data$priors$cointegration$v_i
p_tau_i <- data$priors$cointegration$p_tau_i

sigma_df_prior <- data$priors$sigma$df # Prior degrees of freedom
sigma_scale_prior <- data$priors$sigma$scale # Prior covariance matrix
sigma_df_post <- tt + sigma_df_prior # Posterior degrees of freedom

# Initial values
beta <- matrix(0, k_w, r)
beta[1:r, 1:r] <- diag(1, r)

sigma_i <- diag(1 / .0001, k)

g_i <- sigma_i

iterations <- data$model$iterations # Number of iterations of the Gibbs sampler
burnin <- data$model$burnin # Number of burn-in draws
draws <- iterations + burnin # Total number of draws

# Data containers
draws_alpha <- matrix(NA, k_alpha, iterations)
draws_beta <- matrix(NA, k_beta, iterations)
draws_pi <- matrix(NA, k * k_w, iterations)
draws_gamma <- matrix(NA, k_gamma, iterations)
draws_sigma <- matrix(NA, k^2, iterations)

# Start Gibbs sampler
for (draw in 1:draws) {
  
  # Draw conditional mean parameters
  temp <- post_coint_kls(y = y, beta = beta, w = w, x = x, sigma_i = sigma_i,
                           v_i = v_i, p_tau_i = p_tau_i, g_i = g_i,
                           gamma_mu_prior = a_mu_prior,
                           gamma_v_i_prior = a_v_i_prior)
  alpha <- temp$alpha
  beta <- temp$beta
  Pi <- temp$Pi
  gamma <- temp$Gamma
  
  # Draw variance-covariance matrix
  u <- y - Pi %*% w - matrix(gamma, k) %*% x
  sigma_scale_post <- solve(tcrossprod(u) + v_i * alpha %*% tcrossprod(crossprod(beta, p_tau_i) %*% beta, alpha))
  sigma_i <- matrix(rWishart(1, sigma_df_post, sigma_scale_post)[,, 1], k)
  sigma <- solve(sigma_i)
  
  # Update g_i
  g_i <- sigma_i
  
  # Store draws
  if (draw > burnin) {
    draws_alpha[, draw - burnin] <- alpha
    draws_beta[, draw - burnin] <- beta
    draws_pi[, draw - burnin] <- Pi
    draws_gamma[, draw - burnin] <- gamma
    draws_sigma[, draw - burnin] <- sigma
  }
}

## ----beta---------------------------------------------------------------------
beta <- apply(t(draws_beta) / t(draws_beta)[, 1], 2, mean) # Obtain means for every row
beta <- matrix(beta, k_w) # Transform mean vector into a matrix
beta <- round(beta, 3) # Round values
dimnames(beta) <- list(dimnames(w)[[1]], NULL) # Rename matrix dimensions

beta # Print

## ----bvec-object--------------------------------------------------------------
# Number of non-deterministic coefficients
k_nondet <- (k_x - 4) * k

# Generate bvec object
bvec_est <- bvec(y = data$data$Y,
                 w = data$data$W,
                 x = data$data$X[, 1:6],
                 x_d = data$data$X[, -(1:6)],
                 Pi = draws_pi,
                 Gamma = draws_gamma[1:k_nondet,],
                 C = draws_gamma[(k_nondet + 1):nrow(draws_gamma),],
                 Sigma = draws_sigma)

## -----------------------------------------------------------------------------
summary(bvec_est)

## ---- eval = FALSE------------------------------------------------------------
#  bvec_est <- draw_posterior(data)

## ----thin---------------------------------------------------------------------
bvec_est <- thin_posterior(bvec_est, thin = 5)

## ----vec2var------------------------------------------------------------------
bvar_form <- bvec_to_bvar(bvec_est)

## -----------------------------------------------------------------------------
summary(bvar_form)

## ----forecasts, fig.width=5.5, fig.height=5.5---------------------------------
# Generate deterministc terms for function predict
new_d <- data$data$X[3 + 1:10, c("const", "season.1", "season.2", "season.3")]

# Genrate forecasts
bvar_pred <- predict(bvar_form, n.ahead = 10, new_d = new_d)

# Plot forecasts
plot(bvar_pred)

## ----feir, fig.width=5.5, fig.height=4.5--------------------------------------
FEIR <- irf(bvar_form, impulse = "R", response = "Dp", n.ahead = 20)

plot(FEIR, main = "Forecast Error Impulse Response", xlab = "Period", ylab = "Response")

## ----fevd-oir, fig.width=5.5, fig.height=4.5----------------------------------
bvar_fevd_oir <- fevd(bvar_form, response = "Dp", n.ahead = 20)

plot(bvar_fevd_oir, main = "OIR-based FEVD of inflation")

