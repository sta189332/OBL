#' Compute Optimal Block Length for Non-overlapping, Overlapping, synthetic-overlapping, appended-overlapping and Circular Block Bootstrap
#'
#' More detailed Description
#'
#' @describeIn This package helps to obtain the optimal block length of a time series data
#'
#' @importFrom forecast forecast auto.arima accuracy
#'
#' @param ts univariate time series
#'
#' @param R number of resample
#'
#' @param methods   "optnbb", "optmbb", "optcbb", "opttmbb", "opttcbb"
#'
#' @importFrom forecast auto.arima forecast accuracy
#'
#' @importFrom foreach `%dopar%` foreach
#'
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
#' @importFrom parallel detectCores makeCluster
#'
#' @importFrom future plan multisession
#'
#' @importFrom utils head tail
#'
#' @importFrom stats embed
#'
#' @example
#' set.seed(289805)
#' ts <- arima.sim(n = 10, model = list(ar = 0.8, order = c(1, 0, 0)), sd = 1)
#' blockboot(ts, 1000)
#'
#' @export
blockboot <- function(ts, R, methods = c("optnbb", "optmbb", "optcbb", "opttmbb", "opttcbb")){
  #To ignore the warnings during usage use the first 2 lines
  options(warn = -1)
  options("getSymbols.warning4.0" = FALSE)
  future::plan(future::multisession)
  n_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cores = n_cores)
  nbb <- function(ts, R){
    n <- length(ts)
    lb <- seq(n - 2) + 1
    z <- 1:length(lb)
    `%dopar%` <- foreach::`%dopar%`
    accuracyyy <- foreach::foreach(z = 1:length(lb), .combine = 'rbind', .packages = c('foreach', 'forecast')) %dopar% {
      l <- lb[z]
      blk <- split(ts, ceiling(seq_along(ts) / l)); blk[lengths(blk) == l]
      set.seed(6)
      res <- sample(blk, replace = TRUE, R)
      res.unlist <- unlist(res, use.names = FALSE)

      train <- utils::head(res.unlist, round(length(res.unlist) * 0.20))

      test <- utils::tail(res.unlist, length(res.unlist) - length(train))

      nfuture <- as.numeric(forecast::forecast(train, model = forecast::auto.arima(train), h = length(test))$mean)

      accuracyy <- forecast::accuracy(nfuture, test)
    }

    row.names(accuracyyy) <- lb
    nbb_rmse <- data.frame(lb, accuracyyy[ ,2])
    #return(nbb_rmse[[1]][which.min(nbb_rmse[[2]])])
    colnames(nbb_rmse) <- c("lb", "RMSE")
    with(nbb_rmse, paste(round(min(nbb_rmse$RMSE), 2), nbb_rmse$lb[which.min(nbb_rmse$RMSE)], sep = ' | '))
  }

  # ---------------------------------------------------------------------------
  mbb <- function(ts, R){
    n <- length(ts)
    lb <- seq(n - 2) + 1
    z <- 1:length(lb)
    `%dopar%` <- foreach::`%dopar%`
    accuracyyy <- foreach::foreach(z = 1:length(lb), .combine  =  'rbind', .packages = c('foreach', 'forecast')) %dopar% {
      l <- lb[z]
      b <- n - l + 1
      blk <- split(t(stats::embed(ts, b))[b:1,], 1:b)  # divides the series into overlapping1 blocks
      ######################################################
      set.seed(6)
      res <- sample(blk, replace = T, R)        # resamples the blocks
      res.unlist <- unlist(res, use.names = FALSE)   # unlist the bootstrap series

      train <- head(res.unlist, round(length(res.unlist) * 0.20))

      test <- tail(res.unlist, length(res.unlist) - length(train))

      nfuture <- as.numeric(forecast::forecast(train, model = forecast::auto.arima(train, parallel = TRUE, stepwise = TRUE, num.cores = NULL), h = length(test))$mean)        # makes the `future` object a vector

      accuracyy <- forecast::accuracy(nfuture, test)      # RETURN ACCURACY
    }

    row.names(accuracyyy) <- lb
    mbb_rmse <- data.frame(lb, accuracyyy[,2])
    #return(mbb1_rmse[[1]][which.min(mbb1_rmse[[2]])])
    colnames(mbb_rmse) <- c("lb", "RMSE")
    with(mbb_rmse, paste(round(min(mbb_rmse$RMSE), 2), mbb_rmse$lb[which.min(mbb_rmse$RMSE)], sep = ' | '))
  }

  # ---------------------------------------------------
  cbb <- function(ts, R){
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
      set.seed(6)
      res <- sample(blk, replace = T, R)        # resamples the blocks
      res.unlist <- unlist(res, use.names = FALSE)   # unlist the bootstrap series

      train <- head(res.unlist, round(length(res.unlist) * 0.20))

      test <- tail(res.unlist, length(res.unlist) - length(train))

      nfuture <- as.numeric(forecast::forecast(train, model = forecast::auto.arima(train, parallel = TRUE, stepwise = TRUE, num.cores = NULL), h = length(test))$mean)        # makes the `future` object a vector

      ACCURACYY <- forecast::accuracy(nfuture, test)      # RETURN ACCURACY

      accuracyy <- forecast::accuracy(nfuture, test)      # RETURN ACCURACY
    }

    row.names(accuracyyy) <- lb
    cbb_rmse <- data.frame(lb, accuracyyy[,2])
    #return(cbb_rmse[[1]][min(cbb_rmse[[2]])])
    colnames(cbb_rmse) <- c("lb", "RMSE")
    with(cbb_rmse, paste(round(min(cbb_rmse$RMSE), 2), cbb_rmse$lb[which.min(cbb_rmse$RMSE)], sep = ' | '))
  }

  # ---------------------------------------------------------------

  tmbb <- function(ts, R){
    n <- length(ts)
    lb <- seq(n - 2) + 1
    z <- 1:length(lb)
    `%dopar%` <- foreach::`%dopar%`
    accuracyyy <- foreach::foreach(z = 1:length(lb), .combine  =  'rbind', .packages = c('foreach', 'forecast')) %dopar% {
      #blocks <- function(l, ov, n) {
      #starts <- unique(sort(c(seq(1, n, l), seq(l - ov + 1, n, l))))
      #ends <- pmin(starts + l - 1, n)

      ## truncate starts and ends to the first num elements
      #num <- match(n, ends)
      #head(data.frame(starts, ends), num)
      #}

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
      set.seed(6)
      res <- sample(blk, replace = T, R)        # resamples the blocks
      res.unlist <- unlist(res, use.names = FALSE)   # unlist the bootstrap series

      train <- head(res.unlist, round(length(res.unlist) * 0.20))

      test <- tail(res.unlist, length(res.unlist) - length(train))

      nfuture <- as.numeric(forecast::forecast(train, model = forecast::auto.arima(train, parallel = TRUE, stepwise = TRUE, num.cores = NULL), h = length(test))$mean)        # makes the `future` object a vector

      accuracyy <- forecast::accuracy(nfuture, test)      # RETURN ACCURACY
    }

    row.names(accuracyyy) <- lb
    tmbb_rmse <- data.frame(lb, accuracyyy[,2])
    #return(mbb2_rmse[[1]][which.min(mbb2_rmse[[2]])])
    colnames(tmbb_rmse) <- c("lb", "RMSE")
    with(tmbb_rmse, paste(round(min(tmbb_rmse$RMSE), 2), tmbb_rmse$lb[which.min(tmbb_rmse$RMSE)], sep = ' | '))
  }

  # -------------------------------------------------------------------------

  tcbb <- function(ts, R){
    n <- length(ts)
    lb <- seq(n - 2) + 1
    z <- 1:length(lb)
    `%dopar%` <- foreach::`%dopar%`
    accuracyyy <- foreach::foreach(z  =  1:length(lb), .combine  =  'rbind', .packages = c('foreach', 'forecast')) %dopar% {
      #blocks <- function(l, ov, n) {
      #starts <- unique(sort(c(seq(1, n, l), seq(l - ov + 1, n, l))))
      #ends <- pmin(starts + l - 1, n)

      ## truncate starts and ends to the first num elements
      #num <- match(n, ends)
      #head(data.frame(starts, ends), num)
      #}

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
      set.seed(6)
      res <- sample(blk, replace = T, R)        # resamples the blocks
      res.unlist <- unlist(res, use.names = FALSE)   # unlist the bootstrap series

      train <- head(res.unlist, round(length(res.unlist) * 0.20))

      test <- tail(res.unlist, length(res.unlist) - length(train))

      nfuture <- as.numeric(forecast::forecast(train, model = forecast::auto.arima(train, parallel = TRUE, stepwise = TRUE, num.cores = NULL), h = length(test))$mean)        # makes the `future` object a vector

      accuracyy <- forecast::accuracy(nfuture, test)      # RETURN ACCURACY
    }

    row.names(accuracyyy) <- lb
    tcbb_rmse <- data.frame(lb, accuracyyy[,2])
    #return(mbb3_rmse[[1]][which.min(mbb3_rmse[[2]])])
    colnames(tcbb_rmse) <- c("lb", "RMSE")
    with(tcbb_rmse, paste(round(min(tcbb_rmse$RMSE), 2), tcbb_rmse$lb[which.min(tcbb_rmse$RMSE)], sep = ' | '))
  }



  output <- c()

  if ("optnbb" %in% methods) {
    output <- c(output, nbb = nbb(ts, R))
  }

  if ("optmbb" %in% methods) {
    output <- c(output, mbb = mbb(ts, R))
  }

  if ("optcbb" %in% methods) {
    output <- c(output, cbb = cbb(ts, R))
  }

  if ("opttmbb" %in% methods) {
    output <- c(output, tmbb = tmbb(ts, R))
  }

  if ("opttcbb" %in% methods) {
    output <- c(output, tcbb = tcbb(ts, R))
  }

  return(output)
  doParallel::stopImplicitCluster(cl)
}
