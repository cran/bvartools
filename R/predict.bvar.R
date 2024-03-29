#' Predict Method for Objects of Class bvar
#' 
#' Forecasting a Bayesian VAR object of class \code{"bvar"} with credible bands.
#' 
#' @param object an object of class \code{"bvar"}, usually, a result of a call to
#' \code{\link{bvar}} or \code{\link{bvec_to_bvar}}.
#' @param n.ahead number of steps ahead at which to predict.
#' @param new_x an object of class \code{ts} of new non-deterministic, exogenous variables.
#' The object must have the same frequency as the time series in \code{object[["x"]]} and must contain
#' at least all necessary observations for the predicted period.
#' @param new_d a matrix of new deterministic variables. Must have \code{n.ahead} rows.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the
#' credible intervals. Defaults to 0.95.
#' @param ... additional arguments.
#' 
#' @details For the VAR model
#' \deqn{A_0 y_t = \sum_{i = 1}^{p} A_{i} y_{t-i} + \sum_{i = 0}^{s} B_{i} x_{t-i} + C D_t + u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)} the function produces \code{n.ahead} forecasts.
#' 
#' @return A time-series object of class \code{"bvarprd"}.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' e1 <- window(e1, end = c(1978, 4))
#' 
#' # Generate model data
#' model <- gen_var(e1, p = 0, deterministic = "const",
#'                  iterations = 100, burnin = 10)
#' # Chosen number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' # Generate forecasts
#' bvar_pred <- predict(object, n.ahead = 10, new_d = rep(1, 10))
#' 
#' # Plot forecasts
#' plot(bvar_pred)
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @export
#' @rdname bvar
predict.bvar <- function(object, ..., n.ahead = 10, new_x = NULL, new_d = NULL, ci = .95) {
  
  # Dev specs
  # n.ahead = 10; new_x = NULL; new_d = NULL; ci = .95
  # new_d <- rep(1, 10)
  
  k <- object[["specifications"]][["dims"]][["K"]]
  tt <- nrow(object[["y"]])
  
  if (is.null(object[["A0"]])) {
    struct <- FALSE
  } else {
    struct <- TRUE
  }
  
  A <- NULL
  p <- 0
  tot <- k
  
  tvp <- object[["specifications"]][["tvp"]]
  
  if (!is.null(object[["A"]])) {
    if (tvp[["A"]]) {
      A <- cbind(A, object[["A"]][[tt]])
      p <- ncol(object[["A"]][[1]]) / k^2
    } else {
      A <- cbind(A, object[["A"]])  
      p <- ncol(object[["A"]]) / k^2
    }
    tot <- tot + k * (p - 1)
  }
  
  if (!is.null(object[["B"]])) {
    
    # Get lags
    s <- object[["specifications"]][["lags"]][["s"]]
    
    if (tvp[["B"]]) {
      A <- cbind(A, object[["B"]][[tt]])
      m <- ncol(object[["B"]][[1]]) / k
    } else {
      A <- cbind(A, object[["B"]])
      m <- ncol(object[["B"]]) / k 
    }
    tot <- tot + m
    if (is.null(new_x)) {
      new_x <- matrix(0, n.ahead, m)
    } else {
      
      # Set position of most recent observation, after which forecasts start
      if (s > 0) {
        first_x <- stats::time(object[["y"]])[tt - s + 1]
      }
      
      # Prepare ts column names if only a vector is provided
      if (is.null(dimnames(new_x))) {
        tsp_temp <- stats::tsp(new_x)
        new_x <- stats::ts(as.matrix(new_x), class = c("mts", "ts", "matrix"))
        stats::tsp(new_x) <- tsp_temp
        dimnames(new_x)[[2]] <- "x"
      }
      
      # Obtain lags of exogenous parameters
      if (s > 0) {
        
        # Trim exogenous data to relevant time  
        new_x <- stats::ts(new_x[stats::time(new_x) >= first_x, ], start = first_x, frequency = stats::tsp(new_x)[3]) 
        
        temp <- new_x
        for (i in 1:s) {
          temp <- cbind(temp, stats::lag(new_x, -i))
        }
        new_x <- stats::na.omit(temp)
        dimnames(new_x)[[2]] <- NULL
      }
    }
    
    # Check if new_x is long enough
    if (NROW(new_x) < n.ahead) {
      stop("Longer time series required for argument 'new_x'.")
    }
    
    # Trim exogenous data to forecast horizon
    new_x <- new_x[1:n.ahead,]
  }
  
  if (!is.null(object[["C"]])) {
    if (tvp[["C"]]) {
      A <- cbind(A, object[["C"]][[tt]])
      n <- ncol(object[["C"]][[1]]) / k
    } else {
      A <- cbind(A, object[["C"]])
      n <- ncol(object[["C"]]) / k 
    }
    tot <- tot + n
    if (is.null(new_d)) {
      new_d <- matrix(0, n.ahead, n)
    }
    if (NROW(new_d) != n.ahead) {
      stop("Length of argument 'new_d' must be equal to 'n.ahead'.")
    }
  }
  
  use_a <- !is.null(A)
  
  # Generate matrix used for prediction ----
  pred <- matrix(NA, tot, n.ahead + 1)
  if (p > 0) {
    pos_y <- 1:(k * p)
    pred[pos_y, 1] <- t(object[["y"]][tt:(tt - p + 1), ])
  } else {
    pos_y <- 1:k
    pred[pos_y, 1] <- t(object[["y"]][tt,])
  }
  
  if (!is.null(new_x)) {
    pos_x <- k * p + 1:m
    tt_pos <- stats::time(object[["x"]]) == stats::time(object[["y"]])[tt] # Necessary for BVEC models
    pred[pos_x, ] <- cbind(matrix(object[["x"]][tt_pos, pos_x]), t(new_x))
  }
  
  if (!is.null(object[["C"]])) {
    pos_d_pred <- (tot - n + 1):tot
    if (p == 0) {
      pos_d_object <- (tot - k - n + 1):(tot - k)
    } else {
      pos_d_object <- pos_d_pred
    }
    tt_pos <- stats::time(object[["x"]]) == stats::time(object[["y"]])[tt] # Necessary for BVEC models
    
    # Use latest observation of model data and add new_d
    pred[pos_d_pred, ] <- cbind(matrix(object[["x"]][tt_pos, pos_d_object], length(pos_d_object)), t(new_d))
  }
  
  draws <- NA
  vars <- c("A0", "A", "B", "C", "Sigma")
  for (i in vars) {
    if (is.na(draws)) {
      if (!is.null(object[[i]])) {
        if (object[["specifications"]][["tvp"]][[i]]) {
          draws <- nrow(object[[i]][[1]])
        } else {
          draws <- nrow(object[[i]]) 
        }
      }   
    }
  }
  
  A0_i <- diag(1, k)
  result <- array(NA, dim = c(k, n.ahead, draws))
  for (draw in 1:draws) {
    
    if (struct) {
      if (tvp[["A0"]]) {
        A0_i <- solve(matrix(object[["A0"]][[tt]][draw, ], k))
      } else {
        A0_i <- solve(matrix(object[["A0"]][draw, ], k)) 
      }
    }

    if (tvp[["Sigma"]]) {
      result[,, draw] <- .draw_forecast(draw, k, p, A0_i, use_a, A, object[["Sigma"]][[tt]], pred)[1:k, -1] 
    } else {
      result[,, draw] <- .draw_forecast(draw, k, p, A0_i, use_a, A, object[["Sigma"]], pred)[1:k, -1]
    }
  }
  
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  temp <- apply(result, c(2, 1) , stats::quantile, probs = c(ci_low, .5, ci_high)) 
  result <- c()
  for (i in 1:k) {
    result <- c(result, list(stats::ts(t(temp[,, i]))))
  }
  names(result) <- dimnames(object$y)[[2]]
  
  if (!is.null(attr(object$y, "ts_info"))) {
    ts_info <- attr(object$y, "ts_info")
    object$y <- stats::ts(object$y, start = ts_info[1], frequency = ts_info[3])
    attr(object$y, "ts_info") <- NULL
    
    ts_temp <- stats::ts(0:n.ahead, start = ts_info[2], frequency = ts_info[3])
    ts_temp <- stats::time(ts_temp)[-1]
    for (i in 1:k) {
      stats::tsp(result[[i]]) <- c(ts_temp[1], ts_temp[length(ts_temp)], ts_info[3])
    }
  } else {
    ts_temp <- stats::tsp(object$y)
    for (i in 1:k) {
      result[[i]] <- stats::ts(result[[i]], start = ts_temp[2] + 1 / ts_temp[3], frequency = ts_temp[3])
    }
  }
  
  result <- list("y" = object[["y"]],
                 "fcst" = result)
  
  class(result) <- c("bvarprd", "list")
  return(result)
}
