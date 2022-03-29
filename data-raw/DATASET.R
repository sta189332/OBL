## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

set.seed(289805)
ts <- arima.sim(n = 10, model = list(ar = 0.8, order = c(1, 0, 0)), sd = 1)
usethis::use_data(ts, internal = TRUE, overwrite = TRUE)
