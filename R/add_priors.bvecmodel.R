#' Add Priors for Vector Error Correction Models
#'
#' Adds prior specifications to a list of models, which was produced by
#' function \code{\link{gen_vec}}.
#'
#' @param object a list, usually, the output of a call to \code{\link{gen_vec}}.
#' @param coef a named list of prior specifications for coefficients that do not 
#' determine the cointegration space. For the default specification all prior means are set to zero
#' and the diagonal elements of the inverse prior variance-covariance matrix are set to 1
#' for coefficients corresponding to non-deterministic terms. For deterministic coefficients the prior
#' variances are set to 10 by \code{v_i_det = 0.1}. The variances need to be specified as precisions,
#' i.e. as inverses of the variances. For further specifications such as the Minnesota prior see 'Details'.
#' @param coint a named list of prior specifications for coefficients determining the
#' cointegration space of VEC models. See 'Details'.
#' @param sigma a named list of prior specifications for the error variance-covariance matrix
#' of the models. For the default specification of an inverse Wishart distribution
#' the prior degrees of freedom are set to the number of endogenous variables plus the rank of
#' the cointegration matrix. The prior variance is to 1. See 'Details'.
#' @param ssvs optional; a named list of prior specifications for the SSVS algorithm.
#' Not allowed for TVP models. See 'Details'.
#' @param bvs optional; a named list of prior specifications for the BVS algorithm. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details The arguments of the function require named lists. Possible
#' specifications are described in the following. Note that it is important to specify the
#' priors in the correct list. Otherwise, the provided specification will be disregarded
#' and default values will be used.
#' 
#' Argument \code{coef} contains the following elements
#' \describe{
#'   \item{\code{v_i}}{a numeric specifying the prior precision of the coefficients. Default is 1.}
#'   \item{\code{v_i_det}}{a numeric specifying the prior precision of coefficients that correspond
#'   to deterministic terms. Default is 0.1.}
#'   \item{\code{const}}{a character specifying the prior mean of coefficients, which correspond
#'   to the intercept. If \code{const = "mean"}, the means of the series of endogenous variables
#'   are used as prior means. If \code{const = "first"}, the first values of the series of endogenous
#'   variables are used as prior means.}
#'   \item{\code{minnesota}}{a list of length 4 containing parameters for the calculation of
#'   the Minnesota prior, where the element names must be \code{kappa0}, \code{kappa1}, \code{kappa2} and \code{kappa3}.
#'   For the endogenous variable \eqn{i} the prior variance of the \eqn{l}th lag of regressor \eqn{j} is obtained as
#'   \deqn{ \frac{\kappa_{0}}{l^2} \textrm{ for own lags of endogenous variables,}} 
#'   \deqn{ \frac{\kappa_{0} \kappa_{1}}{l^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for endogenous variables other than own lags,}}
#'   \deqn{ \frac{\kappa_{0} \kappa_{2}}{(l+1)^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for exogenous variables,}}
#'   \deqn{ \kappa_{0} \kappa_{3} \sigma_{i}^2 \textrm{ for deterministic terms,}}
#'   where \eqn{\sigma_{i}} is the residual standard deviation of variable \eqn{i} of an unrestricted
#'   LS estimate. For exogenous variables \eqn{\sigma_{i}} is the sample standard deviation.
#'   
#'   The function only provides priors for the non-cointegration part of the model. However,
#'   the residual standard errors \eqn{\sigma_i} are based on an unrestricted LS regression of the
#'   endogenous variables on the error correction term and the non-cointegration regressors.}
#'   \item{\code{max_var}}{a numeric specifying the maximum prior variance that is allowed for
#'   non-deterministic coefficients.}
#'   \item{\code{shape}}{a numeric specifying the prior shape parameter of the error term of the
#'   state equation. Only used for models with time varying parameters. Default is 3.}
#'   \item{\code{rate}}{a numeric specifying the prior rate parameter of the error term of the
#'   state equation. Only used for models with time varying parameters. Default is 0.0001.}
#'   \item{\code{rate_det}}{a numeric specifying the prior rate parameter of the error term of the
#'   state equation for coefficients, which correspond to deterministic terms.
#'   Only used for models with time varying parameters. Default is 0.01.}
#' }
#' If \code{minnesota} is specified, elements \code{v_i} and \code{v_i_det} are ignored.
#' 
#' Argument \code{coint} can contain the following elements:
#' \describe{
#'   \item{\code{v_i}}{numeric between 0 and 1 specifying the shrinkage of the cointegration space prior. Default is 0.}
#'   \item{\code{p_tau_i}}{numeric of the diagonal elements of the inverse prior matrix of
#'   the central location of the cointegration space \eqn{sp(\beta)}. Default is 1.}
#'   \item{\code{shape}}{an integer specifying the prior degrees of freedom of the error term of the state equation
#'   of the loading matrix \eqn{\alpha}. Default is 3.}
#'   \item{\code{rate}}{a numeric specifying the prior variance of error term of the state equation
#'   of the loading matrix \eqn{\alpha}. Default is 0.0001.}
#'   \item{\code{rho}}{a numeric specifying the autocorrelation coefficient
#'   of the state equation of \eqn{\beta}. It must be smaller than 1. Default is 0.999.
#'   Note that in contrast to Koop et al. (2011) \eqn{\rho} is not drawn in the Gibbs sampler of
#'   this package yet.}
#' }
#' 
#' Argument \code{sigma} can contain the following elements:
#' \describe{
#'   \item{\code{df}}{an integer or character specifying the prior degrees of freedom of the error term. Only
#'   used, if the prior is inverse Wishart.
#'   Default is \code{"k"}, which indicates the amount of endogenous variables in the respective model.
#'   \code{"k + 3"} can be used to set the prior to the amount of endogenous variables plus 3. If an integer
#'   is provided, the degrees of freedom are set to this value in all models.
#'   In all cases the rank \eqn{r} of the cointegration matrix is automatically added.}
#'   \item{\code{scale}}{a numeric specifying the prior error variance of endogenous variables.
#'   Default is 1.}
#'   \item{\code{shape}}{a numeric or character specifying the prior shape parameter of the error term. Only
#'   used, if the prior is inverse gamma or if time varying volatilities are estimated.
#'   For models with constant volatility the default is \code{"k"}, which indicates the amount of endogenous
#'   variables in the respective country model. \code{"k + 3"} can be used to set the prior to the amount of
#'   endogenous variables plus 3. If a numeric is provided, the shape parameters are set to this value in all
#'   models. For models with stochastic volatility this prior refers to the error variance of the state
#'   equation.}
#'   \item{\code{rate}}{a numeric specifying the prior rate parameter of the error term. Only used, if the
#'   prior is inverse gamma or if time varying volatilities are estimated. For models with stochastic
#'   volatility this prior refers to the error variance of the state equation.}
#'   \item{\code{mu}}{numeric of the prior mean of the initial state of the log-volatilities.
#'   Only used for models with time varying volatility.}
#'   \item{\code{v_i}}{numeric of the prior precision of the initial state of the log-volatilities.
#'   Only used for models with time varying volatility.}
#'   \item{\code{sigma_h}}{numeric of the initial draw for the variance of the log-volatilities.
#'   Only used for models with time varying volatility.}
#'   \item{\code{constant}}{numeric of the constant, which is added before taking the log of the squared errors.
#'   Only used for models with time varying volatility.}
#'   \item{\code{covar}}{logical indicating whether error covariances should be estimated. Only used
#'   in combination with an inverse gamma prior or stochastic volatility, for which \code{shape} and
#'   \code{rate} must be specified.}
#' }
#' \code{df} and \code{scale} must be specified for an inverse Wishart prior. \code{shape} and \code{rate}
#' are required for an inverse gamma prior. For structural models or models with stochastic volatility
#' only a gamma prior specification is allowed.
#' 
#' Argument \code{ssvs} can contain the following elements:
#' \describe{
#'   \item{\code{inprior}}{a numeric between 0 and 1 specifying the prior probability
#'   of a variable to be included in the model.}
#'   \item{\code{tau}}{a numeric vector of two elements containing the prior standard errors
#'   of restricted variables (\eqn{\tau_0}) as its first element and unrestricted variables (\eqn{\tau_1})
#'   as its second.}
#'   \item{\code{semiautomatic}}{an numeric vector of two elements containing the
#'   factors by which the standard errors associated with an unconstrained least squares
#'   estimate of the model are multiplied to obtain the prior standard errors
#'   of restricted (\eqn{\tau_0}) and unrestricted (\eqn{\tau_1}) variables, respectively.
#'   This is the semiautomatic approach described in George et al. (2008).}
#'   \item{\code{covar}}{logical indicating if SSVS should also be applied to the error covariance matrix
#'   as in George et al. (2008).}
#'   \item{\code{exclude_det}}{logical indicating if deterministic terms should be excluded from the SSVS algorithm.}
#'   \item{\code{minnesota}}{a numeric vector of length 4 containing parameters for the calculation of
#'   the Minnesota-like inclusion priors. See below.}
#' }
#' Either \code{tau} or \code{semiautomatic} must be specified.
#' 
#' The argument \code{bvs} can contain the following elements
#' \describe{
#'   \item{\code{inprior}}{a numeric between 0 and 1 specifying the prior probability
#'   of a variable to be included in the model.}
#'   \item{\code{covar}}{logical indicating if BVS should also be applied to the error covariance matrix.}
#'   \item{\code{exclude_det}}{logical indicating if deterministic terms should be excluded from the BVS algorithm.}
#'   \item{\code{minnesota}}{a numeric vector of length 4 containing parameters for the calculation of
#'   the Minnesota-like inclusion priors. See below.}
#' }
#' 
#' If either \code{ssvs$minnesota} or \code{bvs$minnesota} is specified, prior inclusion probabilities
#' are calculated in a Minnesota-like fashion as
#' \tabular{cl}{
#' \eqn{\frac{\kappa_1}{l}} \tab for own lags of endogenous variables, \cr
#' \eqn{\frac{\kappa_2}{l}} \tab for other endogenous variables, \cr
#' \eqn{\frac{\kappa_3}{1 + l}} \tab for exogenous variables, \cr
#' \eqn{\kappa_{4}} \tab for deterministic variables, 
#' }
#' for lag \eqn{l} with \eqn{\kappa_1}, \eqn{\kappa_2}, \eqn{\kappa_3}, \eqn{\kappa_4} as the first, second,
#' third and forth element in \code{ssvs$minnesota} or \code{bvs$minnesota}, respectively.
#' 
#' @return A list of models.
#' 
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \doi{10.1016/j.jeconom.2007.08.017}
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
#' simulation for cointegrated models with priors on the cointegration space.
#' \emph{Econometric Reviews, 29}(2), 224--242.
#' \doi{10.1080/07474930903382208}
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2011). Bayesian inference in
#' a time varying cointegration model. \emph{Journal of Econometrics, 165}(2), 210--220.
#' \doi{10.1016/j.jeconom.2011.07.007}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @examples 
#' 
#' # Get data
#' data("e6")
#' 
#' # Create model
#' model <- gen_vec(e6, p = 4, r = 1,
#'                  const = "unrestricted", seasonal = "unrestricted",
#'                  iterations = 100, burnin = 10)
#' # Chosen number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' @export
add_priors.bvecmodel <- function(object,
                                 coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 0.0001, rate_det = 0.01),
                                 coint = list(v_i = 0, p_tau_i = 1, shape = 3, rate = 0.0001, rho = .999),
                                 sigma = list(df = "k", scale = 1, mu = 0, v_i = 0.01, sigma_h = 0.05, constant = 0.0001),
                                 ssvs = NULL,
                                 bvs = NULL,
                                 ...){
  
  only_one_model <- FALSE
  # If only one model is provided, make it compatible with the rest
  if ("data" %in% names(object)) {
    object <- list(object)
    only_one_model <- TRUE
  }
  
  # Checks - Coefficient priors ----
  if (!is.null(coef)) {
    if (!is.null(coef[["v_i"]])) {
      if (coef[["v_i"]] < 0) {
        stop("Argument 'v_i' must be at least 0.")
      } 
      # Define "v_i_det" if not specified (needed for a check later)
      if (is.null(coef[["v_i_det"]])) {
        coef[["v_i_det"]] <- coef[["v_i"]]
      }
    } else {
      if (!any(c("minnesota", "ssvs") %in% names(coef))) {
        stop("If 'coef$v_i' is not specified, at least 'coef$minnesota' or 'coef$ssvs' must be specified.")
      }
    }
  }
  
  if (!is.null(coef[["const"]])) {
    if ("character" %in% class(coef[["const"]])) {
      if (!coef[["const"]] %in% c("first", "mean")) {
        stop("Invalid specificatin of coef$const.")
      }
    }
  }
  
  if (length(coint) < 2) {
    stop("Argument 'coint' must be at least of length 2.")
  }
  
  # Checks - Error priors ----
  if (length(sigma) < 2) {
    stop("Argument 'sigma' must be at least of length 2.")
  } else {
    error_prior <- NULL
    if (any(unlist(lapply(object, function(x) {x[["model"]][["sv"]]})))) { # Check for SV
      if (any(!c("mu", "v_i", "shape", "rate") %in% names(sigma))) {
        stop("Missing prior specifications for stochastic volatility prior.")
      }
      error_prior <- "sv"
    } else {
      if (all(c("shape", "rate") %in% names(sigma))) {
        error_prior <- "gamma"
      }
      if (all(c("df", "scale") %in% names(sigma))) {
        error_prior <- "wishart"
      }
      if (is.null(error_prior)) {
        stop("Invalid specification for argument 'sigma'.")
      }
      if (error_prior == "wishart" & any(unlist(lapply(object, function(x) {x$model$structural})))) {
        stop("Structural models may not use a Wishart prior. Consider specifying argument 'sigma$shape' instead.")
      }
      
      if (error_prior == "wishart") {
        if (sigma[["df"]] < 0) {
          stop("Argument 'sigma$df' must be at least 0.")
        }
        if (sigma[["scale"]] <= 0) {
          stop("Argument 'sigma$scale' must be larger than 0.")
        } 
      }
      if (error_prior == "gamma") {
        if (sigma[["shape"]] < 0) {
          stop("Argument 'sigma$shape' must be at least 0.")
        }
        if (sigma[["rate"]] <= 0) {
          stop("Argument 'sigma$rate' must be larger than 0.")
        } 
      }
    }
  }
  
  # Check Minnesota ----
  minnesota <- FALSE # Minnesota prior?
  if (!is.null(coef[["minnesota"]])) {
    minnesota <- TRUE
  }
  
  # Check SSVS ----
  use_ssvs <- FALSE
  use_ssvs_error <- FALSE
  use_ssvs_semi <- FALSE
  if (!is.null(ssvs)) {
    if (is.null(ssvs[["inprior"]])) {
      stop("Argument 'ssvs$inprior' must be specified for SSVS.")
    }
    if (is.null(ssvs[["tau"]]) & is.null(ssvs[["semiautomatic"]])) {
      stop("Either argument 'ssvs$tau' or 'ssvs$semiautomatic' must be specified for SSVS.")
    }
    if (is.null(ssvs[["exclude_det"]])) {
      ssvs[["exclude_det"]] <- FALSE
    }
    # In case ssvs is specified, check if the semi-automatic approach of 
    # George et al. (2008) should be used
    if (!is.null(ssvs[["semiautomatic"]])) {
      use_ssvs_semi <- TRUE
    }
    
    use_ssvs <- TRUE
    if (minnesota) {
      minnesota <- FALSE
      warning("Minnesota prior specification overwritten by SSVS.")
    }
    
    if (!is.null(ssvs[["covar"]])) {
      if (ssvs[["covar"]]) {
        if (error_prior == "wishart") {
          stop("If SSVS should be applied to error covariances, argument 'sigma$shape'  and 'sigma$rate' must be specified.")
        }
        use_ssvs_error <- TRUE 
      }
      if (is.null(ssvs[["tau"]])) {
        stop("If SSVS should be applied to error covariances, argument 'ssvs$tau' must be specified.")
      }
    }
  }
  
  # Prior a la Korobilis 2013
  use_bvs <- FALSE
  use_bvs_error <- FALSE
  if (!is.null(bvs)) {
    use_bvs <- TRUE
    if (is.null(bvs$inprior)) {
      stop("If BVS should be applied, argument 'bvs$inprior' must be specified.")
    }
    if (is.null(bvs$exclude_det)) {
      bvs$exclude_det <- FALSE
    }
    if (!is.null(bvs$covar)) {
      if (bvs$covar) {
        if (error_prior == "wishart") {
          stop("If BVS should be applied to error covariances, argument 'sigma$shape' must be specified.")
        }
        use_bvs_error <- TRUE 
      }
    }
    if (coef[["v_i"]] == 0 | (coef[["v_i_det"]] == 0 & !bvs[["exclude_det"]])) {
      warning("Using BVS with an uninformative prior is not recommended.")
    }
  }
  
  if (use_ssvs & use_bvs) {
    stop("SSVS and BVS cannot be applied at the same time.")
  }
  
  if (error_prior == "wishart" & (use_ssvs_error | use_bvs_error)) {
    stop("Wishart prior not allowed when BVS or SSVS are applied to covariance matrix.")
  }
  
  varsel_covar <- use_ssvs_error | use_bvs_error
  
  # Generate priors for each country ----
  for (i in 1:length(object)) {
    
    # Get model specs to obtain total number of coeffs
    k <- length(object[[i]][["model"]][["endogen"]][["variables"]])
    p <- object[[i]][["model"]][["endogen"]][["lags"]]
    
    if (k == 1 & (use_ssvs_error | use_bvs_error)) {
      stop("BVS or SSVS cannot be applied to covarianc matrix when there is only one endogenous variable.")
    } 
    
    use_exo <- FALSE
    if (!is.null(object[[i]]$model$exogen)) {
      use_exo <- TRUE
      m <- length(object[[i]]$model$exogen$variables)
      s <- object[[i]]$model$exogen$lags
    } else {
      s <- 0
      m <- 0
    }
    
    # Substract lag from domestic model for VEC
    p <- p - 1
    
    # Total # of non-deterministic coefficients
    n_gamma <- k * (k * p)
    n_upsilon <- k * (m * s)
    
    # Add number of unrestricted deterministic terms
    n_det <- 0
    if (!is.null(object[[i]][["model"]][["deterministic"]])){
      if (!is.null(object[[i]][["model"]][["deterministic"]][["unrestricted"]])) {
        n_det <- length(object[[i]][["model"]][["deterministic"]][["unrestricted"]]) * k
      }
    }
    
    tot_par <- n_gamma + n_upsilon + n_det
    
    covar <- FALSE
    if (!is.null(sigma[["covar"]])) {
      covar <- sigma[["covar"]]
    }
    structural <- object[[i]][["model"]][["structural"]]
    if (covar & structural) {
      stop("Error covariances and structural coefficients cannot be estimated at the same time.")
    }
    sv <- object[[i]][["model"]][["sv"]]
    n_struct <- 0
    n_z <- NCOL(object[[i]][["data"]][["SUR"]])
    if (object[[i]][["model"]][["rank"]] > 0) {
      n_w <- NCOL(object[[i]][["data"]][["W"]])
      n_z <- n_z - k * n_w
    }
    
    if (structural & k > 1) {
      n_struct <- (k - 1) * k / 2
      tot_par <- tot_par + n_struct
    }
    
    #### Cointegration (constant) ----
    n_ect <- k * (k + m)
    if (!is.null(object[[i]][["model"]][["deterministic"]][["restricted"]])) {
      n_ect <- n_ect + k * length(object[[i]][["model"]][["deterministic"]][["restricted"]])
    }
    
    r_temp <- object[[i]][["model"]][["rank"]]
    if (r_temp > 0) {
      
      n_alpha <- r_temp * k
      n_beta <- r_temp * n_ect / k
      
      if (object[[i]][["model"]][["tvp"]]) {
        rho <- coint[["rho"]]
        if (is.null(rho)) {
          rho <- .999
        } else {
          if (rho >= 1) {
            stop("rho must be smaller than 1.")
          }
          if (rho < .8) {
            warning("Value of rho appears very small.")
          }
        }
        object[[i]][["priors"]][["cointegration"]]$alpha <- list(mu = matrix(0, n_alpha),
                                                                      v_i = diag(1 / (1 - rho^2), n_alpha),
                                                                      shape = matrix(coint[["shape"]], n_alpha),
                                                                      rate = matrix(coint[["rate"]], n_alpha))
        object[[i]][["priors"]][["cointegration"]]$beta <- list(mu = matrix(0, n_beta),
                                                                     v_i = diag(1 - rho^2, n_beta))
      } else {
        object[[i]][["priors"]]$cointegration <- list(v_i = coint[["v_i"]],
                                                      p_tau_i = diag(coint[["p_tau_i"]], n_ect / k)) 
      }
    } else {
      if (r_temp == 0 & tot_par == 0) {
        warning("Model with zero cointegration rank and no non-cointegration regressors is skipped.")
      }
    }
    
    #### Non-cointegration (constant) ####
    # Generate prior matrices ----
    if (tot_par > 0) {
      
      # Zero prior means
      mu <- matrix(rep(0, tot_par - n_struct), k)
      
      # Prior for intercept terms
      if (n_det > 0) {
        
        if (!is.null(coef[["const"]]))  {
          
          pos <- which(dimnames(object[[i]]$data$X)[[2]] == "const")
          
          if (length(pos) == 1) {
            if ("character" %in% class(coef[["const"]])) {
              if (coef$const == "first") {
                mu[, pos] <- object[[i]]$data$Y[1, ]
              }
              if (coef$const == "mean") {
                mu[, pos] <- colMeans(object[[i]]$data$Y)
              }
            }
            if ("numeric" %in% class(coef[["const"]])) {
              mu[, pos] <- coef$const
            } 
          }
        }
      }
      
      mu <- matrix(mu)
      
      if (structural) {
        mu <- rbind(mu, matrix(0, n_struct))
      }
      
      object[[i]][["priors"]][["noncointegration"]] <- list(mu = mu)
      
      # Prior covariances
      if (minnesota) {
        # Minnesota prior ----
        minn <- minnesota_prior(object = object[[i]],
                                kappa0 = coef[["minnesota"]][["kappa0"]],
                                kappa1 = coef[["minnesota"]][["kappa1"]],
                                kappa2 = coef[["minnesota"]][["kappa2"]],
                                kappa3 = coef[["minnesota"]][["kappa3"]],
                                max_var = NULL,
                                coint_var = FALSE,
                                sigma = "AR")
        
        object[[i]][["priors"]][["noncointegration"]][["v_i"]] <- minn[["v_i"]]
      }
      
      # SSVS prior ----
      if (use_ssvs) {
        
        if (object[[i]][["model"]][["tvp"]]) {
          stop("SSVS is not supported for TVP models.")
        }
        
        ssvs_temp <- ssvs_prior(object[[i]], tau = ssvs$tau, semiautomatic = ssvs$semiautomatic)
        temp <- inclusion_prior(object[[i]], prob = ssvs$inprior, exclude_deterministics = ssvs$exclude_det,
                                minnesota_like = !is.null(ssvs$minnesota), kappa = ssvs$minnesota)
        object[[i]]$model$varselect <- "SSVS"
        
        # Drop ECT output from inclusion and SSVS output
        if (r_temp > 0) {
          pos_omit <- 1:(k * NCOL(object[[i]]$data$W))
          ssvs_temp$tau0 <- matrix(ssvs_temp$tau0[-pos_omit,])
          ssvs_temp$tau1 <- matrix(ssvs_temp$tau1[-pos_omit,])
          temp$prior <- matrix(temp$prior[-pos_omit,])
          pos_incl <- which(temp$include %in% pos_omit)
          if (length(pos_incl) > 0) {
            temp$include <- matrix(temp$include[-pos_incl])
          }
          temp$include <- temp$include - (temp$include[1] - 1)
          rm(pos_omit)
          rm(pos_incl)
        }
        object[[i]][["priors"]][["noncointegration"]][["v_i"]] <- diag(1 / ssvs_temp$tau1[, 1]^2, tot_par)
        object[[i]][["priors"]][["noncointegration"]][["ssvs"]][["inprior"]] <- temp$prior
        object[[i]][["priors"]][["noncointegration"]][["ssvs"]][["include"]] <- temp$include
        object[[i]][["priors"]][["noncointegration"]][["ssvs"]][["tau0"]] <- ssvs_temp$tau0
        object[[i]][["priors"]][["noncointegration"]][["ssvs"]][["tau1"]] <- ssvs_temp$tau1
        
        rm(temp)
      }
      
      # Regular prior ----
      if (!minnesota & !use_ssvs) {
        v_i <- diag(coef[["v_i"]], tot_par)
        # Add priors for deterministic terms if they were specified
        if (n_det > 0 & !is.null(coef[["v_i_det"]])) {
          diag(v_i)[tot_par - n_struct - n_det + 1:n_det] <- coef[["v_i_det"]]
        }
        
        object[[i]][["priors"]][["noncointegration"]][["v_i"]] <- v_i
      }
      
      if (use_bvs) {
        object[[i]]$model$varselect <- "BVS"
        temp <- inclusion_prior(object[[i]], prob = bvs[["inprior"]], exclude_deterministics = bvs[["exclude_det"]],
                                minnesota_like = !is.null(bvs[["minnesota"]]), kappa = bvs[["minnesota"]])
        
        # Drop ECT output from inclusion and BVS output
        if (r_temp > 0) {
          pos_omit <- 1:(k * NCOL(object[[i]]$data$W)) # Positions of Pi values
          temp[["prior"]] <- matrix(temp[["prior"]][-pos_omit,]) # Incl prior w/o Pi values
          pos_incl <- which(temp[["include"]] %in% pos_omit) # Drop Pi values
          if (length(pos_incl) > 0) {
            temp[["include"]] <- matrix(temp[["include"]][-pos_incl])
          }
          temp[["include"]] <- temp[["include"]] - (temp[["include"]][1] - 1)
          rm(pos_omit)
          rm(pos_incl)  
        }
        object[[i]][["priors"]][["noncointegration"]][["bvs"]][["inprior"]] <- temp[["prior"]]
        object[[i]][["priors"]][["noncointegration"]][["bvs"]][["include"]] <- temp[["include"]]
      }
      
      
      # TVP prior ----
      if (object[[i]][["model"]][["tvp"]]) {
        object[[i]][["priors"]][["noncointegration"]][["shape"]] <- matrix(coef[["shape"]], tot_par)
        object[[i]][["priors"]][["noncointegration"]][["rate"]] <- matrix(coef[["rate"]], tot_par)
        if (n_det > 0 & !is.null(coef[["rate_det"]])) {
          object[[i]][["priors"]][["noncointegration"]][["rate"]][tot_par - n_struct - n_det + 1:n_det, ] <- coef[["rate_det"]]
        }
      }
    }
    
    ## Covar priors ----
    if (!structural & covar & k > 1) {
      
      n_covar <- k * (k - 1) / 2
      object[[i]][["priors"]][["psi"]][["mu"]] <- matrix(0, n_covar)
      object[[i]][["priors"]][["psi"]][["v_i"]] <- diag(coef[["v_i"]], n_covar)
      if (object[[i]][["model"]][["tvp"]]) {
        object[[i]][["priors"]][["psi"]][["shape"]] <- matrix(coef[["shape"]], n_covar)
        object[[i]][["priors"]][["psi"]][["rate"]] <- matrix(coef[["rate"]], n_covar) 
      }
      
      if (use_ssvs_error) {
        object[[i]][["priors"]][["psi"]][["ssvs"]][["inprior"]] <- matrix(ssvs[["inprior"]], n_covar)
        object[[i]][["priors"]][["psi"]][["ssvs"]][["include"]] <- matrix(1:n_covar)
        object[[i]][["priors"]][["psi"]][["ssvs"]][["tau0"]] <- matrix(ssvs[["tau"]][1], n_covar)
        object[[i]][["priors"]][["psi"]][["ssvs"]][["tau1"]] <- matrix(ssvs[["tau"]][2], n_covar)
      }
      
      if (use_bvs_error) {
        object[[i]][["priors"]][["psi"]][["bvs"]][["inprior"]] <- matrix(bvs[["inprior"]], n_covar)
        object[[i]][["priors"]][["psi"]][["bvs"]][["include"]] <- matrix(1:n_covar)
      }
      
    }
    
    # Error term ----
    if (sv) {
      
      object[[i]][["priors"]][["sigma"]][["mu"]] <- matrix(sigma[["mu"]], k)
      object[[i]][["priors"]][["sigma"]][["v_i"]] <- diag(sigma[["v_i"]], k)
      object[[i]][["priors"]][["sigma"]][["shape"]] <- matrix(sigma[["shape"]], k)
      object[[i]][["priors"]][["sigma"]][["rate"]] <- matrix(sigma[["rate"]], k)
      
    } else {
      if (error_prior == "wishart") {
        object[[i]][["priors"]][["sigma"]][["type"]] <- "wishart"
        help_df <- sigma[["df"]]
        object[[i]][["priors"]][["sigma"]]$df <- NA
        object[[i]][["priors"]][["sigma"]]$scale = diag(sigma[["scale"]], k)
      }
      if (error_prior == "gamma") {
        object[[i]][["priors"]][["sigma"]][["type"]] <- "gamma"
        help_df <- sigma[["shape"]]
        object[[i]][["priors"]][["sigma"]]$shape <- NA
        object[[i]][["priors"]][["sigma"]]$rate = matrix(sigma[["rate"]], k)
      }
      if (minnesota) {
        # Store LS estimate of variance coviariance matrix for analytical solution
        object[[i]][["priors"]][["sigma"]]$sigma_i = minn[["sigma_i"]]
      }
      
      if ("character" %in% class(help_df)) {
        if (grepl("k", help_df)) {
          # Transform character specification to expression and evaluate
          help_df <- eval(parse(text = help_df))
        } else {
          stop("Use no other letter than 'k' in 'sigma$df' to indicate the number of endogenous variables.")
        }
      }
      
      if (help_df < 0) {
        stop("Current specification implies a negative prior degree of\nfreedom or shape parameter of the error term.")
      }
      
      # Add rank to degrees of freedom for cointegration model
      if (!is.na(object[[i]]$model$rank)) {
        help_df <- help_df + r_temp
      }
      
      if (error_prior == "wishart") {
        object[[i]][["priors"]][["sigma"]]$df <- help_df
      }
      if (error_prior == "gamma") {
        object[[i]][["priors"]][["sigma"]]$shape <- matrix(help_df, k)
      }
    }
    
    # Initial values ----
    y <- t(object[[i]][["data"]][["Y"]])
    if (tot_par < NCOL(y)) {
      z <- object[[i]][["data"]][["SUR"]]
      ols <- solve(crossprod(z)) %*% crossprod(z, matrix(y))
      u <- matrix(matrix(y) - z %*% ols, NROW(y))
      object[[i]][["initial"]][["noncointegration"]]$draw <- matrix(ols[-c(1:n_ect),])
    } else {
      if (tot_par > 0) {
        object[[i]][["initial"]][["noncointegration"]]$draw <- mu 
      }
      u <- y - matrix(apply(y, 1, mean), nrow = NROW(y), ncol = NCOL(y))
    }
    if (r_temp > 0) {
      beta <- matrix(0, n_ect / k, r_temp)
      beta[1:r_temp, 1:r_temp] <- diag(1, r_temp)
      object[[i]][["initial"]][["cointegration"]]$beta <- beta
      if (object[[i]][["model"]][["tvp"]]) {
        object[[i]][["initial"]][["cointegration"]][["rho"]] <- rho
      }
    }
    if (object[[i]][["model"]][["tvp"]]) {
      object[[i]][["initial"]][["noncointegration"]][["sigma_i"]] <- diag(c(1 / object[[i]][["priors"]][["noncointegration"]][["rate"]]), tot_par)
    }
    if (covar) {
      y_covar <- kronecker(-t(u), diag(1, k))
      pos <- NULL
      for (j in 1:k) {pos <- c(pos, (j - 1) * k + 1:j)}
      y_covar <- y_covar[, -pos]
      psi <- solve(crossprod(y_covar)) %*% crossprod(y_covar, matrix(u))
      object[[i]][["initial"]][["psi"]][["draw"]] <- psi
      Psi <- diag(1, k)
      for (j in 2:k) {
        Psi[j, 1:(j - 1)] <- t(psi[((j - 2) * (j - 1) / 2) + 1:(j - 1), 1])
      }
      u <- Psi %*% u
    }
    u <- apply(u, 1, stats::var)
    if (object[[i]][["model"]][["sv"]]) {
      object[[i]][["initial"]][["sigma"]][["h"]] <- log(matrix(u, nrow = NCOL(y), ncol = NROW(y), byrow = TRUE))
      object[[i]][["initial"]][["sigma"]][["sigma_h"]] <- matrix(sigma[["sigma_h"]], NROW(y))
      
      if (is.null(sigma[["constant"]])) {
        warning("Argument 'sigma$constant' not specified. Using the value 0.0001.")
        sigma[["constant"]] <- .0001
      }
      object[[i]][["initial"]][["sigma"]][["constant"]] <- matrix(sigma[["constant"]], NROW(y))
    } else {
      object[[i]][["initial"]][["sigma"]][["sigma_i"]] <- diag(1 / u, NROW(y))
      dimnames(object[[i]][["initial"]][["sigma"]][["sigma_i"]]) <- list(dimnames(y)[[2]], dimnames(y)[[2]])
    }
    
  } # End model loop
  
  if (only_one_model) {
    object <- object[[1]]
  }
  
  return(object)
}