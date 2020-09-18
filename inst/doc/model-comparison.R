## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE, warning = FALSE, fig.height=4, fig.width=5-------------
library(bvartools)

data("e1") # Load data
data <- diff(log(e1)) * 100 # Obtain log-differences

# Use date up to 1978Q4
data <- window(data, end = c(1978, 4))

# Plot
plot(data)

## -----------------------------------------------------------------------------
object <- gen_var(data, p = 0:4,
                  deterministic = "const",
                  iterations = 5000, burnin = 1000)

## -----------------------------------------------------------------------------
object <- add_priors(object,
                     coef = list(v_i = 0, v_i_det = 0),
                     sigma = list(df = "k", scale = 0.0001))

## ---- message = FALSE, warning=FALSE, results='hide', eval = FALSE------------
#  object <- draw_posterior(object, mc.cores = 3)

## ---- message = FALSE, warning=FALSE, results='hide', echo = FALSE------------
object <- draw_posterior(object)

## -----------------------------------------------------------------------------
summary(object)

## -----------------------------------------------------------------------------
plot(irf(object[[3]], impulse = "income", response = "cons", n.ahead = 10))

