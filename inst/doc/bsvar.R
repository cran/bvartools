## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data, fig.align='center', fig.height=5, fig.width=4.5---------------
library(bvartools)

data("e1")
e1 <- diff(log(e1))

plot(e1) # Plot the series

## ------------------------------------------------------------------------
data <- gen_var(e1, p = 2, deterministic = "const")

y <- data$Y[, 1:73]
x <- data$Z[, 1:73]

## ----flat prior----------------------------------------------------------
# Reset random number generator for reproducibility
set.seed(1234567)

iter <- 10000 # Number of iterations of the Gibbs sampler
burnin <- 5000 # Number of burn-in draws
store <- iter - burnin

t <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
m <- k * nrow(x) # Number of estimated coefficients
k0 <- k * (k - 1) / 2 # Number of structural coefficients

# Set (uninformative) priors
a_mu_prior <- matrix(0, m) # Vector of prior parameter means
a_v_i_prior <- diag(0, m) # Inverse of the prior covariance matrix

a0_mu_prior <- matrix(0, k0) # Vector of prior parameter means
a0_v_i_prior <- diag(0, k0) # Inverse of the prior covariance matrix

sigma_df_prior <- 0 # Prior degrees of freedom
sigma_scale_prior <- rep(0, k)
sigma_df_post <- t + sigma_df_prior

# Initial values
a0 <- diag(1, k)
omega_i <- rWishart(1, t, solve(tcrossprod(y)))[,, 1]
omega <- solve(omega_i)
sigma_i <- diag(1, k)
diag(sigma_i) <- diag(omega_i)
sigma <- solve(sigma_i)

# Data containers for posterior draws
draws_a <- matrix(NA, m, store)
draws_a0 <- matrix(NA, k^2, store)
draws_omega <- matrix(NA, k^2, store)

# Start Gibbs sampler
for (draw in 1:iter) {
  # Draw conditional mean parameters
  a <- post_normal(y, x, omega_i, a_mu_prior, a_v_i_prior)
  
  # Structural coefficients
  y_tilde <- y - matrix(a, k) %*% x # Obtain residuals
  for (j in 2:k) {
   # Preparing the data
   y_tilde_temp <- matrix(y_tilde[j, ], 1)
   x0_temp <- matrix(-y_tilde[1:(j - 1),], j - 1)
   a0_sigma_i_temp <-  matrix(sigma_i[j, j])
   pos_temp <- (j - 1) * (j - 2) / 2 + 1:(j - 1)
   mu_temp <- matrix(a0_mu_prior[pos_temp,])
   v_i_temp <- matrix(a0_v_i_prior[pos_temp, pos_temp], j - 1)
   
   # Draw structural coefficients
   a0_temp <- post_normal(y = y_tilde_temp, x = x0_temp, sigma_i = a0_sigma_i_temp,
                          a_prior = mu_temp, v_i_prior = v_i_temp)
   
   # Update A0 matrix
   a0[j, 1:(j - 1)] <- a0_temp
  }
  
  # Draw variances
  y_star <- a0 %*% y_tilde
  sigma_scale_post <- sigma_scale_prior + rowSums(y_star^2)
  for (j in 1:k) {
    sigma_i[j, j] <- rgamma(1, shape = sigma_df_post / 2, rate = sigma_scale_post[j] / 2)
  }
  sigma <- solve(sigma_i)
  
  a0_i <- solve(a0)
  omega <- a0_i %*% tcrossprod(sigma, a0_i) 
  omega_i <- solve(omega)
  
  # Store draws
  if (draw > burnin) {
    draws_a[, draw - burnin] <- a
    draws_a0[, draw - burnin] <- a0
    draws_omega[, draw - burnin] <- sigma
  }
}

## ------------------------------------------------------------------------
A <- rowMeans(draws_a) # Obtain means for every row
A <- matrix(A, k) # Transform mean vector into a matrix
A <- round(A, 3) # Round values
dimnames(A) <- list(dimnames(y)[[1]], dimnames(x)[[1]]) # Rename matrix dimensions

A # Print

## ------------------------------------------------------------------------
A0 <- rowMeans(draws_a0) # Obtain means for every row
A0 <- matrix(A0, k) # Transform mean vector into a matrix
A0 <- round(A0, 3) # Round values
dimnames(A0) <- list(dimnames(y)[[1]], dimnames(y)[[1]]) # Rename matrix dimensions

solve(A0) # Print

## ----bvar-object---------------------------------------------------------
bvar_est <- bvar(y = y, x = x, A = draws_a[1:18,], C = draws_a[19:21, ],
                 A0 = draws_a0, Sigma = draws_omega)

## ----thin----------------------------------------------------------------
bvar_est <- thin(bvar_est, thin = 5)

## ----forecasts, fig.width=5.5, fig.height=5.5----------------------------
bvar_pred <- predict(bvar_est, n.ahead = 10, new_D = rep(1, 10))

plot(bvar_pred)

## ----oir, fig.width=5.5, fig.height=4.5----------------------------------
SIR <- irf(bvar_est, impulse = "income", response = "cons", n.ahead = 8, type = "sir")

plot(SIR, main = "Structural Impulse Response", xlab = "Period", ylab = "Response")

## ----fevd-oir, fig.width=5.5, fig.height=4.5-----------------------------
SFEVD <- fevd(bvar_est, response = "cons", n.ahead = 8, type = "sir")

plot(SFEVD, main = "Structrual FEVD of consumption")

