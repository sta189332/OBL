#' @title OBL: Optimal Block Length
#'
#' Compute Optimal Block Length for Non-overlapping, Overlapping, Circular Block, tapered moving, and tapered circular Block Bootstrap method
#'
#' @describeIn This package helps to obtain the optimal block length of a time series data
#'
#' @importFrom forecast forecast auto.arima accuracy
#'
#' @param ts univariate time series
#'
#' @param R number of resample
#'
#' @param seed RNG seed
#'
#' @param n_cores number of core(s) to be used on your operaterating system
#'
#' @param methods   "optnbb", "optmbb", "optcbb", "opttmbb", "opttcbb"
#'
#' @importFrom forecast auto.arima forecast accuracy
#'
#' @importFrom forcats fct_reorder
#'
#' @importFrom foreach `%dopar%` foreach
#'
#' @importFrom dplyr mutate
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_segment scale_color_manual theme_bw
#'
#' @importFrom utils head tail
#'
#' @importFrom stats embed
#'
#' @importFrom tibble rownames_to_column
#'
#' @return A data frame get printed to the console
#'
#' @examples
#'   set.seed(289805)
#'   ts <- arima.sim(n = 3, model = list(ar = 0.8, order = c(1, 0, 0)), sd = 1)
#'   lolliblock(ts, R = 2, seed = 6, n_cores = 1)
#'
#' @export
lolliblock <- function(ts, R, seed, n_cores, methods = c("optnbb", "optmbb", "optcbb", "opttmbb", "opttcbb")){
  #To ignore the warnings during usage use the first 2 lines
  suppressWarnings('non')#options(warn = -1)
  options("getSymbols.warning4.0" = FALSE)
  #future::plan(future::multisession)
  #nn_cores <- parallel::detectCores()
  #cl <- parallel::makeCluster(nn_cores)
  ##on.exit(parallel::stopCluster(cl))
  #doParallel::registerDoParallel(cores = n_cores)
  nbb <- function(ts, R, seed, n_cores){
    n <- length(ts)
    lb <- seq(n - 2) + 1
    z <- 1:length(lb)
    `%dopar%` <- foreach::`%dopar%`
    accuracyyy <- foreach::foreach(z = 1:length(lb), .combine = 'rbind', .packages = c('foreach', 'forecast')) %dopar% {
      l <- lb[z]
      blk <- split(ts, ceiling(seq_along(ts) / l)); blk[lengths(blk) == l]
      set.seed(seed)
      res <- sample(blk, replace = TRUE, R)
      res.unlist <- unlist(res, use.names = FALSE)

      train <- utils::head(res.unlist, round(length(res.unlist) - 2))

      test <- utils::tail(res.unlist, length(res.unlist) - length(train))

      nfuture <- as.numeric(forecast::forecast(train, model = forecast::auto.arima(train, parallel = TRUE, stepwise = FALSE, num.cores = n_cores), h = length(test))$mean)

      accuracyy <- forecast::accuracy(nfuture, test)
    }

    row.names(accuracyyy) <- lb
    nbb_rmse <- data.frame(lb, accuracyyy[ ,2])
    colnames(nbb_rmse) <- c("lb", "RMSE")
    nbb_rmse
  }

  # ---------------------------------------------------------------------------
  mbb <- function(ts, R, seed, n_cores){
    n <- length(ts)
    lb <- seq(n - 2) + 1
    z <- 1:length(lb)
    `%dopar%` <- foreach::`%dopar%`
    accuracyyy <- foreach::foreach(z = 1:length(lb), .combine  =  'rbind', .packages = c('foreach', 'forecast')) %dopar% {
      l <- lb[z]
      b <- n - l + 1
      blk <- split(t(stats::embed(ts, b))[b:1,], 1:b)  # divides the series into overlapping1 blocks
      ######################################################
      set.seed(seed)
      res <- sample(blk, replace = T, R)        # resamples the blocks
      res.unlist <- unlist(res, use.names = FALSE)   # unlist the bootstrap series

      train <- head(res.unlist, round(length(res.unlist) - 2))

      test <- tail(res.unlist, length(res.unlist) - length(train))

      nfuture <- as.numeric(forecast::forecast(train, model = forecast::auto.arima(train, parallel = TRUE, stepwise = FALSE, num.cores = n_cores), h = length(test))$mean)        # makes the `future` object a vector

      accuracyy <- forecast::accuracy(nfuture, test)      # RETURN ACCURACY
    }

    row.names(accuracyyy) <- lb
    mbb_rmse <- data.frame(lb, accuracyyy[,2])
    colnames(mbb_rmse) <- c("lb", "RMSE")
    mbb_rmse
  }

  # ---------------------------------------------------
  cbb <- function(ts, R, seed, n_cores){
    n <- length(ts)
    lb <- seq(n - 2) + 1
    z <- 1:length(lb)
    `%dopar%` <- foreach::`%dopar%`
    accuracyyy <- foreach::foreach(z = 1:length(lb), .combine  =  'rbind', .packages = c('foreach', 'forecast')) %dopar% {
      l <- lb[z]
      newts <- c(ts, head(ts, n = 1))
      m <- length(newts) - l + 1
      blk <- split(t(stats::embed(newts, m))[m:1,], 1:m)
      ######################################################
      set.seed(seed)
      res <- sample(blk, replace = T, R)        # resamples the blocks
      res.unlist <- unlist(res, use.names = FALSE)   # unlist the bootstrap series

      train <- head(res.unlist, round(length(res.unlist) - 2))

      test <- tail(res.unlist, length(res.unlist) - length(train))

      nfuture <- as.numeric(forecast::forecast(train, model = forecast::auto.arima(train, parallel = TRUE, stepwise = FALSE, num.cores = n_cores), h = length(test))$mean)        # makes the `future` object a vector

      ACCURACYY <- forecast::accuracy(nfuture, test)      # RETURN ACCURACY

      accuracyy <- forecast::accuracy(nfuture, test)      # RETURN ACCURACY
    }

    row.names(accuracyyy) <- lb
    cbb_rmse <- data.frame(lb, accuracyyy[,2])
    colnames(cbb_rmse) <- c("lb", "RMSE")
    cbb_rmse
  }

  # ---------------------------------------------------------------

  tmbb <- function(ts, R, seed, n_cores){
    n <- length(ts)
    lb <- seq(n - 2) + 1
    z <- 1:length(lb)
    `%dopar%` <- foreach::`%dopar%`
    accuracyyy <- foreach::foreach(z = 1:length(lb), .combine  =  'rbind', .packages = c('foreach', 'forecast')) %dopar% {

      blockss <- function(l, ov, n) {
        starts <- pmin(n - l + 1, unique(sort(c(seq(1, n, l), seq(l - ov + 1, n, l)))))
        starts <- unique(sort(c(seq(1, n, l), seq(l - ov + 1, n, l))))
        ends <- pmin(starts + l - 1, n)

        # truncate starts and ends to the first num elements
        num <- match(n, ends)
        head(data.frame(starts, ends), num)
      }
      l <- lb[z]                                          # block size
      ov = ceiling(l/2)
      d <- blockss(l, ov, n)
      blk <- with(d, Map(function(i, j) ts[i:j], starts, ends))  # divides the series into overlapping1 blocks
      ######################################################
      set.seed(seed)
      res <- sample(blk, replace = T, R)        # resamples the blocks
      res.unlist <- unlist(res, use.names = FALSE)   # unlist the bootstrap series

      train <- head(res.unlist, round(length(res.unlist) - 2))

      test <- tail(res.unlist, length(res.unlist) - length(train))

      nfuture <- as.numeric(forecast::forecast(train, model = forecast::auto.arima(train, parallel = TRUE, stepwise = FALSE, num.cores = n_cores), h = length(test))$mean)        # makes the `future` object a vector

      accuracyy <- forecast::accuracy(nfuture, test)      # RETURN ACCURACY
    }

    row.names(accuracyyy) <- lb
    tmbb_rmse <- data.frame(lb, accuracyyy[,2])
    #return(mbb2_rmse[[1]][which.min(mbb2_rmse[[2]])])
    colnames(tmbb_rmse) <- c("lb", "RMSE")
    tmbb_rmse
  }

  # -------------------------------------------------------------------------

  tcbb <- function(ts, R, seed, n_cores){
    n <- length(ts)
    lb <- seq(n - 2) + 1
    z <- 1:length(lb)
    `%dopar%` <- foreach::`%dopar%`
    accuracyyy <- foreach::foreach(z  =  1:length(lb), .combine  =  'rbind', .packages = c('foreach', 'forecast')) %dopar% {
      blockss <- function(l, ov, n) {
        starts <- pmin(n - l + 1, unique(sort(c(seq(1, n, l), seq(l - ov + 1, n, l)))))
        starts <- unique(sort(c(seq(1, n, l), seq(l - ov + 1, n, l))))
        ends <- pmin(starts + l - 1, n)

        # truncate starts and ends to the first num elements
        num <- match(n, ends)
        head(data.frame(starts, ends), num)
      }
      l <- lb[z]                                          # block size
      ov = ceiling(l/2)
      d <- blockss(l, ov, n)
      #m <- ceiling(n / l)
      blk <- Map(function(i) {x <- seq(i, i + l - 1) %% max(ts + 1); x + cumsum(x == 0)}, d$starts) # or
      ######################################################
      set.seed(seed)
      res <- sample(blk, replace = T, R)        # resamples the blocks
      res.unlist <- unlist(res, use.names = FALSE)   # unlist the bootstrap series

      train <- head(res.unlist, round(length(res.unlist) - 2))

      test <- tail(res.unlist, length(res.unlist) - length(train))

      nfuture <- as.numeric(forecast::forecast(train, model = forecast::auto.arima(train, parallel = TRUE, stepwise = FALSE, num.cores = n_cores), h = length(test))$mean)        # makes the `future` object a vector

      accuracyy <- forecast::accuracy(nfuture, test)      # RETURN ACCURACY
    }
    row.names(accuracyyy) <- lb
    tcbb_rmse <- data.frame(lb, accuracyyy[,2])
    colnames(tcbb_rmse) <- c("lb", "RMSE")
    tcbb_rmse
  }

  output <- list()

  if ("optnbb" %in% methods) {
    output <- c(output, nbb = nbb(ts, R, seed, n_cores))
  }

  if ("optmbb" %in% methods) {
    output <- c(output, mbb = mbb(ts, R, seed, n_cores))
  }

  if ("optcbb" %in% methods) {
    output <- c(output, cbb = cbb(ts, R, seed, n_cores))
  }

  if ("opttmbb" %in% methods) {
    output <- c(output, tmbb = tmbb(ts, R, seed, n_cores))
  }

  if ("opttcbb" %in% methods) {
    output <- c(output, tcbb = tcbb(ts, R, seed, n_cores))
  }

  df <- list(nbb = data.frame(output$nbb.lb, output$nbb.RMSE), mbb = data.frame(output$mbb.lb, output$mbb.RMSE), cbb = data.frame(output$cbb.lb, output$cbb.RMSE), tmbb = data.frame(output$tmbb.lb, output$tmbb.RMSE), tcbb = data.frame(output$tcbb.lb, output$tcbb.RMSE))
  df1 <- do.call(rbind, lapply(df, function(x) data.frame(lb = x[which.min(x[,2]), 1], RMSE = min(x[, 2])))) |>
    tibble::rownames_to_column("Methods")

  df1 |>
    dplyr::mutate(colour = forcats::fct_reorder(df1$Methods, df1$RMSE)) |>
    ggplot2::ggplot(ggplot2::aes(df1$Methods, df1$RMSE, colour = df1$colour)) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_segment(ggplot2::aes(df1$Methods, xend = df1$Methods, yend = df1$RMSE, y = 0)) +
    ggplot2::scale_color_manual(values = c("green", "yellowgreen", "yellow",
                                           "orange", "red"),
                                labels = c(9, 8, 9, 9, 4), name = "lb") +
    ggplot2::theme_bw(base_size = 16)
}
