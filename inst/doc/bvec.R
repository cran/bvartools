## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data, fig.align='center', fig.height=5, fig.width=4.5---------------
library(bvartools)

data("e6")

plot(e6) # Plot the series

## ------------------------------------------------------------------------
data <- gen_vec(e6, p = 4, const = "unrestricted", season = "unrestricted")

y <- data$Y
w <- data$W
x <- data$X

## ----flat prior----------------------------------------------------------
# Reset random number generator for reproducibility
set.seed(1234567)

iter <- 10000 # Number of iterations of the Gibbs sampler
burnin <- 5000 # Number of burn-in draws
store <- iter - burnin

r <- 1 # Set rank

t <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
k_w <- nrow(w) # Number of regressors in error correction term
k_x <- nrow(x) # Number of differenced regressors and unrestrictec deterministic terms

k_alpha <- k * r # Number of elements in alpha
k_beta <- k_w * r # Number of elements in beta
k_gamma <- k * k_x

# Set uninformative priors
a_mu_prior <- matrix(0, k_x * k) # Vector of prior parameter means
a_v_i_prior <- diag(0, k_x * k) # Inverse of the prior covariance matrix

v_i <- 0
p_tau_i <- diag(1, k_w)

u_sigma_df_prior <- r # Prior degrees of freedom
u_sigma_scale_prior <- diag(0, k) # Prior covariance matrix
u_sigma_df_post <- t + u_sigma_df_prior # Posterior degrees of freedom

# Initial values
beta <- matrix(c(1, -4), k_w, r)
#beta <- matrix(rnorm(k_w * r), k_w, r)

u_sigma_i <- diag(.0001, k)
u_sigma <- solve(u_sigma_i)

g_i <- u_sigma_i

# Data containers
draws_alpha <- matrix(NA, k_alpha, store)
draws_beta <- matrix(NA, k_beta, store)
draws_pi <- matrix(NA, k * k_w, store)
draws_gamma <- matrix(NA, k_gamma, store)
draws_sigma <- matrix(NA, k^2, store)

# Start Gibbs sampler
for (draw in 1:iter) {
  # Draw conditional mean parameters
  temp <- post_coint_kls(y = y, beta = beta, w = w, x = x, sigma_i = u_sigma_i,
                           v_i = v_i, p_tau_i = p_tau_i, g_i = g_i,
                           gamma_mu_prior = a_mu_prior,
                           gamma_V_i_prior = a_v_i_prior)
  alpha <- temp$alpha
  beta <- temp$beta
  Pi <- temp$Pi
  gamma <- temp$Gamma
  
  # Draw variance-covariance matrix
  u <- y - Pi %*% w - matrix(gamma, k) %*% x
  u_sigma_scale_post <- solve(tcrossprod(u) + v_i * alpha %*% tcrossprod(crossprod(beta, p_tau_i) %*% beta, alpha))
  u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
  u_sigma <- solve(u_sigma_i)
  
  # Update g_i
  g_i <- u_sigma_i
  
  # Store draws
  if (draw > burnin) {
    draws_alpha[, draw - burnin] <- alpha
    draws_beta[, draw - burnin] <- beta
    draws_pi[, draw - burnin] <- Pi
    draws_gamma[, draw - burnin] <- gamma
    draws_sigma[, draw - burnin] <- u_sigma
  }
}

## ------------------------------------------------------------------------
Gamma <- rowMeans(draws_gamma) # Obtain means for every row
Gamma <- matrix(Gamma, k) # Transform mean vector into a matrix
Gamma <- round(Gamma, 3) # Round values
dimnames(Gamma) <- list(dimnames(y)[[1]], dimnames(x)[[1]]) # Rename matrix dimensions

round(Gamma, 3) # Print

## ----beta----------------------------------------------------------------
beta <- rowMeans(t(t(draws_beta) / t(draws_beta)[, 1])) # Obtain means for every row
beta <- matrix(beta, k_w) # Transform mean vector into a matrix
beta <- round(beta, 3) # Round values
dimnames(beta) <- list(dimnames(w)[[1]], NULL) # Rename matrix dimensions

beta # Print

## ------------------------------------------------------------------------
Sigma <- rowMeans(draws_sigma) # Obtain means for every row
Sigma <- matrix(Sigma, k) # Transform mean vector into a matrix
Sigma <- round(Sigma * 10^4, 2) # Round values
dimnames(Sigma) <- list(dimnames(y)[[1]], dimnames(y)[[1]]) # Rename matrix dimensions

Sigma # Print

## ----bvec-object---------------------------------------------------------
# Number of non-deterministic coefficients
k_nondet <- (k_x - 4) * k

# Generate bvec object
bvec_est <- bvec(y = y, w = w, x = x,
               Pi = draws_pi,
               Gamma = draws_gamma[1:k_nondet,],
               C = draws_gamma[(k_nondet + 1):nrow(draws_gamma),],
               Sigma = draws_sigma)

## ----thin----------------------------------------------------------------
bvec_est <- thin(bvec_est, thin = 5)

## ----vec2var-------------------------------------------------------------
bvar_form <- bvec_to_bvar(bvec_est)

## ----forecast, , fig.width=5.5, fig.height=5.5---------------------------
bvar_pred <- predict(bvar_form, n.ahead = 10, new_D = t(bvar_form$x[9:12, 91:100]))

plot(bvar_pred)

## ----feir, fig.width=5.5, fig.height=4.5---------------------------------
IR <- irf(bvar_form, impulse = "R", response = "Dp", n.ahead = 20)

plot(IR, main = "Forecast Error Impulse Response", xlab = "Year", ylab = "Response")

## ----oir, fig.width=5.5, fig.height=4.5----------------------------------
OIR <- irf(bvar_form, impulse = "R", response = "Dp", n.ahead = 20, type = "oir")

plot(OIR, main = "Orthogonalised Impulse Response", xlab = "Year", ylab = "Response")

## ----gir, fig.width=5.5, fig.height=4.5----------------------------------
GIR <- irf(bvar_form, impulse = "R", response = "Dp", n.ahead = 20, type = "gir")

plot(GIR, main = "Generalised Impulse Response", xlab = "Year", ylab = "Response")

## ----fevd, fig.width=5.5, fig.height=4.5---------------------------------
bvec_fevd <- fevd(bvar_form, response = "Dp", n.ahead = 20)

plot(bvec_fevd, main = "FEVD of inflation")
