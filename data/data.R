set.seed(289805)
ts <- stats::arima.sim(n = 10, model = list(ar = 0.8, order = c(1, 0, 0)), sd = 1)
