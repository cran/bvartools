#' Bayesian Vector Autoregression Objects
#' 
#' \code{bvar} is used to create objects of class \code{"bvar"}.
#' 
#' @param data the original time-series object of endogenous variables.
#' @param exogen the original time-series object of unmodelled variables.
#' @param y a time-series object of endogenous variables with \eqn{T} observations,
#' usually, a result of a call to \code{\link{gen_var}}.
#' @param x a time-series object of \eqn{(pK + (1+s)M + N)} regressor variables, usually, a result of a
#' call to \code{\link{gen_var}}.
#' @param A0 either a \eqn{K^2 \times S} matrix of MCMC coefficient draws of structural parameters or
#' a named list, where element \code{coeffs} contains a \eqn{K^2 \times S} matrix of MCMC coefficient
#' draws of structural parameters and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed. For time varying parameter models
#' the coefficient matrix must be \eqn{TK^2 \times S}. Draws of the error covariance matrix of the state
#' equation can be provided as a \eqn{K^2 \times S} matrix in an additional list element.
#' @param A either a \eqn{pK^2 \times S} matrix of MCMC coefficient draws of lagged endogenous variables or
#' a named list, where element \code{coeffs} contains a \eqn{pK^2 \times S} matrix of MCMC coefficient draws
#' of lagged endogenous variables and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed. For time varying parameter models
#' the coefficient matrix must be \eqn{pTK^2 \times S}. Draws of the error covariance matrix of the state
#' equation can be provided as a \eqn{pK^2 \times S} matrix in an additional list element.
#' @param B either a \eqn{((1 + s)MK) \times S} matrix of MCMC coefficient draws of unmodelled, non-deterministic variables
#' or a named list, where element \code{coeffs} contains a \eqn{((1 + s)MK) \times S} matrix of MCMC coefficient draws of
#' unmodelled, non-deterministic variables and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed. For time varying parameter models
#' the coefficient matrix must be \eqn{(1 + s)TMK \times S}. Draws of the error covariance matrix of the state
#' equation can be provided as a \eqn{(1 + s)MK \times S} matrix in an additional list element.
#' @param C either a \eqn{KN \times S} matrix of MCMC coefficient draws of deterministic terms or
#' a named list, where element \code{coeffs} contains a \eqn{KN \times S} matrix of MCMC coefficient draws of
#' deterministic terms and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed. For time varying parameter models
#' the coefficient matrix must be \eqn{TKN \times S}. Draws of the error covariance matrix of the state
#' equation can be provided as a \eqn{KN \times S} matrix in an additional list element.
#' @param Sigma a \eqn{K^2 \times S} matrix of MCMC draws for the error variance-covariance matrix or
#' a named list, where element \code{coeffs} contains a \eqn{K^2 \times S} matrix of MCMC draws for the
#' error variance-covariance matrix and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed to the covariances. For time varying parameter models
#' the coefficient matrix must be \eqn{TK^2 \times S}. Draws of the error covariance matrix of the state
#' equation can be provided as a \eqn{K^2 \times S} matrix in an additional list element.
#' 
#' @details For the VARX model
#' \deqn{A_0 y_t = \sum_{i = 1}^{p} A_i y_{t-i} + \sum_{i = 0}^{s} B_i x_{t - i} + C d_t + u_t}
#' the function collects the S draws of a Gibbs sampler (after the burn-in phase) in a standardised object,
#' where \eqn{y_t} is a K-dimensional vector of endogenous variables,
#' \eqn{A_0} is a \eqn{K \times K} matrix of structural coefficients.
#' \eqn{A_i} is a \eqn{K \times K} coefficient matrix of lagged endogenous variabels.
#' \eqn{x_t} is an M-dimensional vector of unmodelled, non-deterministic variables
#' and \eqn{B_i} its corresponding coefficient matrix.
#' \eqn{d_t} is an N-dimensional vector of deterministic terms
#' and \eqn{C} its corresponding coefficient matrix.
#' \eqn{u_t} is an error term with \eqn{u_t \sim N(0, \Sigma_u)}.
#' 
#' For time varying parameter and stochastic volatility models the respective coefficients and
#' error covariance matrix of the above model are assumed to be time varying, respectively.
#' 
#' The draws of the different coefficient matrices provided in \code{A0}, \code{A},
#' \code{B}, \code{C} and \code{Sigma} have to correspond to the same MCMC iterations.
#' 
#' @return An object of class \code{"bvar"} containing the following components, if specified:
#' \item{data}{the original time-series object of endogenous variables.}
#' \item{exogen}{the original time-series object of unmodelled variables.}
#' \item{y}{a \eqn{K \times T} matrix of endogenous variables.}
#' \item{x}{a \eqn{(pK + (1+s)M + N) \times T} matrix of regressor variables.}
#' \item{A0}{an \eqn{S \times K^2} "mcmc" object of coefficient draws of structural parameters. In case of time varying parameters a list of such objects.}
#' \item{A0_lambda}{an \eqn{S \times K^2} "mcmc" object of inclusion parameters for structural parameters.}
#' \item{A0_sigma}{an \eqn{S \times K^2} "mcmc" object of the error covariance matrices of the structural parameters in a model with time varying parameters.}
#' \item{A}{an \eqn{S \times pK^2} "mcmc" object of coefficient draws of lagged endogenous variables. In case of time varying parameters a list of such objects.}
#' \item{A_lambda}{an \eqn{S \times pK^2} "mcmc" object of inclusion parameters for lagged endogenous variables.}
#' \item{A_sigma}{an \eqn{S \times pK^2} "mcmc" object of the error covariance matrices of coefficients of lagged endogenous variables in a model with time varying parameters.}
#' \item{B}{an \eqn{S \times ((1 + s)MK)} "mcmc" object of coefficient draws of unmodelled, non-deterministic variables. In case of time varying parameters a list of such objects.}
#' \item{B_lambda}{an \eqn{S \times ((1 + s)MK)} "mcmc" object of inclusion parameters for unmodelled, non-deterministic variables.}
#' \item{B_sigma}{an \eqn{S \times ((1 + s)MK)} "mcmc" object of the error covariance matrices of coefficients of unmodelled, non-deterministic variables in a model with time varying parameters.}
#' \item{C}{an \eqn{S \times NK} "mcmc" object of coefficient draws of deterministic terms. In case of time varying parameters a list of such objects.}
#' \item{C_lambda}{an \eqn{S \times NK} "mcmc" object of inclusion parameters for deterministic terms.}
#' \item{C_sigma}{an \eqn{S \times NK} "mcmc" object of the error covariance matrices of coefficients of deterministic terms in a model with time varying parameters.}
#' \item{Sigma}{an \eqn{S \times K^2} "mcmc" object of variance-covariance draws. In case of time varying parameters a list of such objects.}
#' \item{Sigma_lambda}{an \eqn{S \times K^2} "mcmc" object of inclusion parameters for error covariances.}
#' \item{Sigma_sigma}{an \eqn{S \times K^2} "mcmc" object of the error covariance matrices of the coefficients of the error covariance matrix of the measurement equation of a model with time varying parameters.}
#' \item{specifications}{a list containing information on the model specification.}

#' @examples
#' 
#' # Get data
#' data("e1")
#' e1 <- diff(log(e1))
#' e1 <- window(e1, end = c(1978, 4))
#' 
#' # Generate model data
#' data <- gen_var(e1, p = 2, deterministic = "const")
#'
#' # Add priors
#' model <- add_priors(data,
#'                     coef = list(v_i = 0, v_i_det = 0),
#'                     sigma = list(df = 0, scale = .00001))
#' 
#' # Set RNG seed for reproducibility 
#' set.seed(1234567)
#' 
#' iterations <- 400 # Number of iterations of the Gibbs sampler
#' # Chosen number of iterations and burnin should be much higher.
#' burnin <- 100 # Number of burn-in draws
#' draws <- iterations + burnin # Total number of MCMC draws
#'
#' y <- t(model$data$Y)
#' x <- t(model$data$Z)
#' tt <- ncol(y) # Number of observations
#' k <- nrow(y) # Number of endogenous variables
#' m <- k * nrow(x) # Number of estimated coefficients
#' 
#' # Priors
#' a_mu_prior <- model$priors$coefficients$mu # Vector of prior parameter means
#' a_v_i_prior <- model$priors$coefficients$v_i # Inverse of the prior covariance matrix
#' 
#' u_sigma_df_prior <- model$priors$sigma$df # Prior degrees of freedom
#' u_sigma_scale_prior <- model$priors$sigma$scale # Prior covariance matrix
#' u_sigma_df_post <- tt + u_sigma_df_prior # Posterior degrees of freedom
#'
#' # Initial values
#' u_sigma_i <- diag(1 / .00001, k)
#'
#' # Data containers for posterior draws
#' draws_a <- matrix(NA, m, iterations)
#' draws_sigma <- matrix(NA, k^2, iterations)
#'
#' # Start Gibbs sampler
#' for (draw in 1:draws) {
#'  # Draw conditional mean parameters
#'  a <- post_normal(y, x, u_sigma_i, a_mu_prior, a_v_i_prior)
#'
#'  # Draw variance-covariance matrix
#'  u <- y - matrix(a, k) %*% x # Obtain residuals
#'  u_sigma_scale_post <- solve(u_sigma_scale_prior + tcrossprod(u))
#'  u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
#'
#'  # Store draws
#'  if (draw > burnin) {
#'   draws_a[, draw - burnin] <- a
#'   draws_sigma[, draw - burnin] <- solve(u_sigma_i)
#'  }
#' }
#' 
#' # Generate bvar object
#' bvar_est <- bvar(y = model$data$Y, x = model$data$Z,
#'                  A = draws_a[1:18,], C = draws_a[19:21, ],
#'                  Sigma = draws_sigma)
#'                  
#' @export
bvar <- function(data = NULL, exogen = NULL, y, x = NULL,
                 A0 = NULL, A = NULL, B = NULL,
                 C = NULL, Sigma = NULL) {
  
  if (!"ts" %in% class(y)) {
    stop("Argument 'y' must be an object of class time-series")
  }
  
  result <- NULL
  result[["y"]] <- y
  k <- NCOL(y)
  tt <- NROW(y)
  
  tvp_a0 <- FALSE
  tvp_a <- FALSE
  tvp_b <- FALSE
  tvp_c <- FALSE
  tvp_sigma <- FALSE
  
  structural <- FALSE
  
  if(!is.null(A0)) {
    if (is.list(A0)) {
      if ("coeffs" %in% names(A0)) {
        n_a0 <- nrow(A0[["coeffs"]])
      }
    } else {
      n_a0 <- nrow(A0)
    }
    if (n_a0 / tt >= 1) {
      tvp_a0 <- TRUE
      n_a0 <- n_a0 / tt
    }
    if (n_a0 %% (k * k) != 0) {
      stop("Row number of coefficient draws of 'A0' is not k^2 or multiples thereof.")
    }
    structural <- TRUE
  }
  
  if(!is.null(A)) {
    if (is.list(A)) {
      if ("coeffs" %in% names(A)) {
        n_a <- nrow(A[["coeffs"]])
      }
    } else {
      n_a <- nrow(A)
    }
    if ((n_a / tt) %% k^2 == 0 & n_a / tt >= 1) {
      tvp_a <- TRUE
      n_a <- n_a / tt
    }
  }
  
  if(!is.null(B)) {
    if (is.null(x) & is.null(exogen)) {
      stop("Please specify either argument 'x' or 'exogen' when using exogenous regressors.")
    }
    if (is.list(B)) {
      if ("coeffs" %in% names(B)) {
        n_b <- nrow(B[["coeffs"]])
      }
    } else {
      n_b <- nrow(B)
    }
    if ((n_b / tt) %% k == 0) {
      tvp_b <- TRUE
      n_b <- n_b / tt
    }
  }
  
  if(!is.null(C)) {
    if (is.list(C)) {
      if ("coeffs" %in% names(C)) {
        n_c <- NROW(C[["coeffs"]])
      }
    } else {
      n_c <- NROW(C)
    }
    if ((n_c / tt) %% k == 0 & n_c / tt >= 1) {
      tvp_c <- TRUE
      n_c <- n_c / tt
    }
  }
  
  if(!is.null(Sigma)) {
    if (is.list(Sigma)) {
      if ("coeffs" %in% names(Sigma)) {
        n_sigma <- nrow(Sigma[["coeffs"]])
      }
    } else {
      n_sigma <- nrow(Sigma)
    }
    if ((n_sigma / tt) %% k == 0 & n_sigma / tt >= 1) {
      tvp_sigma <- TRUE
      n_sigma <- n_sigma / tt
    }
    if (n_sigma %% (k * k) != 0) {
      stop("Row number of coefficient draws of 'Sigma' is not k^2 or multiples thereof.")
    }
  }
  
  # Data objects ----
  if(!is.null(data)) {
    result[["data"]] <- data
  }
  if(!is.null(exogen)) {
    result[["exogen"]] <- exogen
  }
  if(!is.null(x)) {
    result[["x"]] <- x
  }
  
  # Parameters - A0 ----
  if(!is.null(A0)) {
    result <- c(result, .bvar_fill_helper(A0, tvp_a0, n_a0, tt, "A0"))
  }
  
  # Parameters - A ----
  if(!is.null(A)) {
    if (n_a %% k == 0) {
      p <- n_a / k^2 
    } else {
      stop("Row number of argument 'A' is not a multiple of the number of endogenous variables.")
    }
    result <- c(result, .bvar_fill_helper(A, tvp_a, n_a, tt, "A"))
  } else {
    p <- 0
  }
  
  # Parameters - B ----
  m <- 0
  s <- 0
  if(!is.null(B)) {
    result <- c(result, .bvar_fill_helper(B, tvp_b, n_b, tt, "B"))
    if (!is.null(exogen)) {
      m <- NCOL(exogen)
      s <- n_b / (k * m) - 1 
    }
    if (!is.null(x)) {
      x_names <- dimnames(x)[[2]][k * p + 1:(n_b / k)]
      x_lags <- as.numeric(substring(x_names, unlist(gregexpr(".l", x_names)) + 2, nchar(x_names)))
      m <- sum(x_lags == 0)
      s <- max(x_lags)
    }
  }
  
  # Parameters - C ----
  if(!is.null(C)) {
    result <- c(result, .bvar_fill_helper(C, tvp_c, n_c, tt, "C"))
  }
  
  # Parameters - Sigma ----
  if(!is.null(Sigma)) {
    result <- c(result, .bvar_fill_helper(Sigma, tvp_sigma, k * k, tt, "Sigma"))
  }
  
  result[["specifications"]] <- list("dims" = list("K" = k, "M" = m),
                                     "lags" = list("p" = p, "s" = s),
                                     "tvp" = list("A0" = tvp_a0,
                                                  "A" = tvp_a,
                                                  "B" = tvp_b,
                                                  "C" = tvp_c,
                                                  "Sigma" = tvp_sigma),
                                     "structural" = structural)
  
  class(result) <- "bvar"
  return(result)
}