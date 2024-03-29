% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/draw_posterior.bvarmodel.R
\name{draw_posterior.bvarmodel}
\alias{draw_posterior.bvarmodel}
\title{Posterior Simulation}
\usage{
\method{draw_posterior}{bvarmodel}(object, FUN = NULL, mc.cores = NULL, ...)
}
\arguments{
\item{object}{a list of model specifications, which should be passed on
to function \code{FUN}. Usually, the output of a call to \code{\link{gen_var}} in
combination with \code{\link{add_priors}}.}

\item{FUN}{the function to be applied to each model in argument \code{object}.
If \code{NULL} (default), the internal functions \code{\link{bvarpost}} is used.}

\item{mc.cores}{the number of cores to use, i.e. at most how many child
processes will be run simultaneously. The option is initialized from
environment variable MC_CORES if set. Must be at least one, and
parallelization requires at least two cores.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
For multiple models a list of objects of class \code{bvarlist}.
For a single model the object has the class of the output of the applied posterior
simulation function. In case the package's own functions are used, this will
result in an object of class \code{"bvar"}.
}
\description{
Forwards model input to posterior simulation functions.
}
\examples{

# Load data 
data("e1")
e1 <- diff(log(e1)) * 100

# Generate model
model <- gen_var(e1, p = 1:2, deterministic = 2,
                 iterations = 100, burnin = 10)
# Chosen number of iterations and burn-in should be much higher.

# Add priors
model <- add_priors(model)

# Obtain posterior draws
object <- draw_posterior(model)

}
