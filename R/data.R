#' Ten (10) simulated univaariate time series data.
#'
#' \code{arima.sim} returns the sum of all the values present in its arguments.
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
#' @source Simulated data generated with the following code:
#'   set.seed(289805)
#'   ts <- stats::arima.sim(n = 10, model = list(ar = 0.8, order = c(1, 0, 0)), sd = 1)
#'
#' @return It returns a univairate time series data
#' It could be a vector
#'
#' @examples  set.seed(289805)
#' @examples  ts <- stats::arima.sim(n = 10, model = list(ar = 0.8, order = c(1, 0, 0)), sd = 1)
"ts"
