#' Ten (10) simulated univaariate time series data.
#'
#' A dataset containing simulated univariate time series of 10
#' ts.
#'
#' @format A time series data with 10 rows and 1 variables:
#' \describe{
#'   \item{price}{price, in US dollars}
#'   \item{carat}{weight of the diamond, in carats}
#'   ...
#' }
#' @source \url{(set.seed(289805)ts <- stats::arima.sim(n = 10, model = list(ar = 0.8, order = c(1, 0, 0)), sd = 1)}
"ts"
