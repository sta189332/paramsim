#' @title Parameterized Simulation
#'
#' @describeIn This arimasim helps to Search for rigth seeds for the rigth AR simulation with arima.sin() finction using auto.arima() function
#'
#' Search for rigth seeds for the rigth ARIMA simulation with arima.sin() function using auto.arima() function
#'
#' This function obtains a Random Number Generator (RNG) or collection of RNGs that replicate the required parameter(s) of a distribution for a time series of data. Consider the case of reproducing a time series data set of size 20 that uses an autoregressive (AR) model with phi = 0.8 and standard deviation equal to 1. When one checks the arima.sin() function's estimated parameters, it's possible that after a single trial or a few more, one won't find the precise parameters. This enables one to look for the ideal RNG setting for a simulation that will accurately duplicate the desired parameters.
#'
#' @param a first seed boundary
#'
#' @param z last seed boundary
#'
#' @param n number of samples
#'
#' @param ar11 first coefficient of autoregressive
#'
#' @param ar22 second coefficient of autoregressive
#'
#' @param ar33 third coefficient of autoregressive
#'
#' @param p order of the autoregressive
#'
#' @param d degree of difference
#'
#' @param  q degree of moving average
#'
#' @param sd standard deviation of the series
#'
#' @param j1 length of character to search for in first coefficient of autoregressive
#'
#' @param j2 length of character to search for in second coefficient of autoregressive
#'
#' @param j3 length of character to search for in third coefficient of autoregressive
#'
#' @param k1 length of character to search for in third coefficient of autoregressive
#'
#' @param k2 length of character to search for in third coefficient of autoregressive
#'
#' @param k3 length of character to search for in third coefficient of autoregressive
#'
#' @param arr1 character to search for in first coefficient of autoregressive
#'
#' @param arr2 character to search for in second coefficient of autoregressive
#'
#' @param arr3 character to search for in third coefficient of autoregressive
#'
#' @param ar11 character to search for in third coefficient of autoregressive
#'
#' @param ar22 character to search for in third coefficient of autoregressive
#'
#' @param ar33 character to search for in third coefficient of autoregressive
#'
#' @param ma11 character to search for in third coefficient of autoregressive
#'
#' @param ma22 character to search for in third coefficient of autoregressive
#'
#' @param ma33 character to search for in third coefficient of autoregressive
#'
#' @param maa1 character to search for in third coefficient of autoregressive
#'
#' @param maa2 character to search for in third coefficient of autoregressive
#'
#' @param maa3 character to search for in third coefficient of autoregressive
#'
#' @importFrom future plan multisession
#'
#' @importFrom parallel makeCluster stopCluster
#'
#' @importFrom doParallel registerDoParallel
#'
#' @importFrom foreach `%dopar%` foreach
#'
#' @importFrom stats arima.sim
#'
#' @importFrom forecast auto.arima
#'
#' @importFrom tibble tibble
#'
#' @return A data frame get printed to the console with its first colomn being the rank and the next few column could be the coefficients of AR or MA both with varying orders depending on the order and classes of ARIMA model being searched for. The last column of the data frame could be the intercept if any exist within the range of the search.
#'
#' @examples
#'   arimasim(a= 289805,z= 289806,n= 10,p= 1,d= 0,q= 0,ar11= 0.8,sd = 1,j1= 4,arr1= "0.80")
#'
#' @export
arimasim <- function(a, z, n, ar11, ma11, ar22, ma22, ar33, ma33, p, d, q, sd = sd, j1, k1, j2, k2, j3, k3, arr1, maa1, arr2, maa2, arr3, maa3){

  output <- if (p == 1 && q == 0) {

    ar1_sim_search <- function(a, z, n, ar11, p, d, q, sd = sd, j1, arr1){

      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ar = c(ar11), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        } else if (all(grepl(c("ar1|intercept"), names(cf))) &
                   substr(cf["ar1"], 1, j1) %in% arr1) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]
      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ar1_sim_search(a = a,  z = z, n = n, p = 1, d = d, q = q, ar11 = ar11, sd = sd, j1 = j1, arr1 = arr1)

    #####################################################################################################

  } else if (p == 0 && q == 1) {

    ma1_sim_search <- function(a, z, n, ma11, p, d, q, sd = sd, k1, maa1){

      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ma = c(ma11), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        } else if (all(grepl(c("ma1|intercept"), names(cf))) &
                   substr(cf["ma1"], 1, k1) %in% maa1) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]
      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ma1_sim_search(a = a,  z = z, n = n, p = p, d = d, q = 1, ma11 = ma11, sd = sd, k1 = k1, maa1 = maa1)

    #####################################################################################################

    #####################################################################################################

  } else if (p == 1 && q == 1) {

    arma1_sim_search <- function(a, z, n, ar11, ma11, p, d, q, sd = sd, j1, k1, arr1, maa1){
      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ar = ar11, ma = ma11, order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        }  else if (all(grepl(c("ar1|ma1|intercept"), names(cf))) &
                    substr(cf["ar1"], 1, j1) %in% arr1 & substr(cf["ma1"], 1, k1) %in% maa1) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]


      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    arma1_sim_search(a = a,  z = z, n = n, p = 1, d = d, q = 1, ar11 = ar11, ma11 = ma11, sd = sd, j1 = j1, k1 = k1, arr1 = arr1, maa1 = maa1)

    #####################################################################################################

  } else if (p == 2 && q == 0) {

    ar2_sim_search <- function(a, z, n, ar11, ar22, p, d, q, sd = sd, j1, j2, arr1, arr2){
      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ar = c(ar11, ar22), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        }  else if (all(grepl(c("ar1|ar2|intercept"), names(cf))) &
                    substr(cf["ar1"], 1, j1) %in% arr1 & substr(cf["ar2"], 1, j2) %in% arr2) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]


      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ar2_sim_search(a = a,  z = z, n = n, p = 2, d = d, q = q, ar11 = ar11, ar22 = ar22, sd = sd, j1 = j1, j2 = j2, arr1 = arr1, arr2 = arr2)

    #####################################################################################################

  } else if (p == 0 && q == 2) {

    ma2_sim_search <- function(a, z, n, ma11, ma22, p, d, q, sd = sd, k1, k2, maa1, maa2){
      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ma = c(ma11, ma22), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        }  else if (all(grepl(c("ma1|ma2|intercept"), names(cf))) &
                    substr(cf["ma1"], 1, k1) %in% maa1 & substr(cf["ma2"], 1, k2) %in% maa2) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]


      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ma2_sim_search(a = a,  z = z, n = n, p = 0, d = 0, q = 2, ma11 = ma11, ma22 = ma22, sd = sd, k1 = k1, k2 = k2, maa1 = maa1, maa2 = maa2)
    ##############################################################################

  } else if (p == 3 && q == 0) {

    ar3_sim_search <- function(a, z, n, ar11, ar22, ar33, p, d, q, sd = sd, j1, j2, j3, arr1, arr2, arr3){
      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ar = c(ar11, ar22, ar33), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        } else if (all(grepl(c("ar1|ar2|ar3|intercept"), names(cf))) &
                   substr(cf["ar1"], 1, j1) %in% arr1 & substr(cf["ar2"], 1, j2) %in% arr2 & substr(cf["ar3"], 1, j3) %in% arr3) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]


      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ar3_sim_search(a = a,  z = z, n = n, p = 3, d = d, q = q, ar11 = ar11, ar22 = ar22, ar33 = ar33, sd = sd, j1 = j1, j2 = j2, j3 = j3, arr1 = arr1, arr2 = arr2, arr3 = arr3)

    #####################################################################################################

  } else if (p == 3 && d == 1 && q == 0) {

    ma3_sim_search <- function(a, z, n, ma11, ma22, ma33, p, d, q, sd = sd, k1, k2, k3, maa1, maa2, maa3){
      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ma = c(ma11, ma22, ma33), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        } else if (all(grepl(c("ma1|ma2|ma3|intercept"), names(cf))) &
                   substr(cf["ma1"], 1, k1) %in% maa1 & substr(cf["ma2"], 1, k2) %in% maa2 & substr(cf["ma3"], 1, k3) %in% maa3) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]


      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ma3_sim_search(a = a,  z = z, n = n, p = p, d = d, q = 3, ma11 = ma11, ma22 = ma22, ma33 = ma33, sd = sd, k1 = k1, k2 = k2, k3 = k3, maa1 = maa1, maa2 = maa2, maa3 = maa3)

    #}
    ##############################################################################
    ##############################################################################
    ##############################################################################
  } else if (p == 1 && d == 1 && q == 0) {

    ari1_sim_search <- function(a, z, n, ar11, p, d, q, sd = sd, j1, arr1){
      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ar = c(ar11), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        } else if (all(grepl(c("ar1|intercept"), names(cf))) &
                   substr(cf["ar1"], 1, j1) %in% arr1) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]


      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ari1_sim_search(a = a,  z = z, n = n, p = 1, d = 1, q = q, ar11 = ar11, sd = sd, j1 = j1, arr1 = arr1)

    #####################################################################################################

  } else if (p == 0 && d == 1 && q == 1) {

    ima1_sim_search <- function(a, z, n, ma11, p, d, q, sd = sd, k1, maa1){
      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ma = c(ma11), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        } else if (all(grepl(c("ma1|intercept"), names(cf))) &
                   substr(cf["ma1"], 1, k1) %in% maa1) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]


      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ima1_sim_search(a = a,  z = z, n = n, p = p, d = 1, q = 1, ma11 = ma11, sd = sd, k1 = k1, maa1 = maa1)

    #####################################################################################################

    #####################################################################################################

  } else if (p == 1 && d == 1 && q == 1) {

    arima11_sim_search <- function(a, z, n, ar11, ma11, p, d, q, sd = sd, j1, k1, arr1, maa1){
      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ar = ar11, ma = ma11, order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        }  else if (all(grepl(c("ar1|ma1|intercept"), names(cf))) &
                    substr(cf["ar1"], 1, j1) %in% arr1 & substr(cf["ma1"], 1, k1) %in% maa1) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]


      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    arima11_sim_search(a = a,  z = z, n = n, p = 1, d = 1, q = 1, ar11 = ar11, ma11 = ma11, sd = sd, j1 = j1, k1 = k1, arr1 = arr1, maa1 = maa1)

    #####################################################################################################

  } else if (p == 2 && d == 1 && q == 0) {

    ari2_sim_search <- function(a, z, n, ar11, ar22, p, d, q, sd = sd, j1, j2, arr1, arr2){
      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ar = c(ar11, ar22), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        }  else if (all(grepl(c("ar1|ar2|intercept"), names(cf))) &
                    substr(cf["ar1"], 1, j1) %in% arr1 & substr(cf["ar2"], 1, j2) %in% arr2) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]


      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ari2_sim_search(a = a,  z = z, n = n, p = 2, d = 1, q = 0, ar11 = ar11, ar22 = ar22, sd = sd, j1 = j1, j2 = j2, arr1 = arr1, arr2 = arr2)

    #####################################################################################################

  } else if (p == 0 && d == 1 && q == 2) {

    ima2_sim_search <- function(a, z, n, ma11, ma22, p, d, q, sd = sd, k1, k2, maa1, maa2){
      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ma = c(ma11, ma22), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        }  else if (all(grepl(c("ma1|ma2|intercept"), names(cf))) &
                    substr(cf["ma1"], 1, k1) %in% maa1 & substr(cf["ma2"], 1, k2) %in% maa2) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]


      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ima2_sim_search(a = a,  z = z, n = n, p = 0, d = 1, q = 2, ma11 = ma11, ma22 = ma22, sd = sd, k1 = k1, k2 = k2, maa1 = maa1, maa2 = maa2)
    ##############################################################################

  } else if (p == 1 && d == 1 && q == 1) {

    arima22_sim_search <- function(a, z, n, ar11, ar22, ma11, ma22, p, d, q, sd = sd, j1, j2, k1, k2, arr1, arr2, maa1, maa2){
      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ar = ar11, ma = ma11, order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        }  else if (all(grepl(c("ar1|ar2|ma1|ma2|intercept"), names(cf))) &
                    substr(cf["ar1"], 1, j1) %in% arr1 & substr(cf["ar2"], 1, j2) %in% arr2 & substr(cf["ma1"], 1, k1) %in% maa1 & substr(cf["ma2"], 1, k2) %in% maa2) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]


      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    arima22_sim_search(a = a,  z = z, n = n, p = 2, d = 1, q = 2, ar11 = ar11, ma11 = ma11, ar22 = ar22, ma22 = ma22, sd = sd, j1 = j1, j2 = j2, k1 = k1, k2 = k2, arr1 = arr1, arr2 = arr2, maa1 = maa1, maa2 = maa2)
    ##############################################################################


  } else if (p == 3 && d == 1 && q == 0) {

    ari3_sim_search <- function(a, z, n, ar11, ar22, ar33, p, d, q, sd = sd, j1, j2, j3, arr1, arr2, arr3){
      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ar = c(ar11, ar22, ar33), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        } else if (all(grepl(c("ar1|ar2|ar3|intercept"), names(cf))) &
                   substr(cf["ar1"], 1, j1) %in% arr1 & substr(cf["ar2"], 1, j2) %in% arr2 & substr(cf["ar3"], 1, j3) %in% arr3) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]


      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ari3_sim_search(a = a,  z = z, n = n, p = 3, d = 1, q = 0, ar11 = ar11, ar22 = ar22, ar33 = ar33, sd = sd, j1 = j1, j2 = j2, j3 = j3, arr1 = arr1, arr2 = arr2, arr3 = arr3)

    #####################################################################################################

  } else if (p == 0 && d == 1 && q == 3) {

    ima3_sim_search <- function(a, z, n, ma11, ma22, ma33, p, d, q, sd = sd, k1, k2, k3, maa1, maa2, maa3){
      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ma = c(ma11, ma22, ma33), order = c(p, d, q)), sd = sd)
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        } else if (all(grepl(c("ma1|ma2|ma3|intercept"), names(cf))) &
                   substr(cf["ma1"], 1, k1) %in% maa1 & substr(cf["ma2"], 1, k2) %in% maa2 & substr(cf["ma3"], 1, k3) %in% maa3) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]


      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    ima3_sim_search(a = a,  z = z, n = n, p = 0, d = 1, q = 3, ma11 = ma11, ma22 = ma22, ma33 = ma33, sd = sd, k1 = k1, k2 = k2, k3 = k3, maa1 = maa1, maa2 = maa2, maa3 = maa3)

    ##############################################################################

  } else {

    arima33_sim_search <- function(a, z, n, ar11, ar22, ar33, ma11, ma22, ma33, p, d, q, sd = 1, j1, j2, j3, k1, k2, k3, arr1, arr2, arr3, maa1, maa2, maa3){
      message('processing...')
      `%dopar%` <- foreach::`%dopar%`
      i <- a:z
      res <- foreach::foreach(i = a:z, .packages = c('foreach', 'forecast')) %dopar% {
        set.seed(i)
        mod <- stats::arima.sim(n = n, model = list(ar = c(ar11, ar22, ar33), ma = c(ma11, ma22, ma33), order = c(p, d, q)), sd = sd)
        #mod <- arima.sim(model = list(ar = c(ar11), ma = c(ma11)), n = n, order = c(p, d, q), n.start = 1, start.innov = rnorm(n, sd = 1), rand.gen = (1 + 2 * ar11 * ma11 + ma11^2) * 1 / (1 - ar11^2))
        best.mod <- forecast::auto.arima(mod, ic = "aicc")
        (cf <- best.mod$coef)
        if (length(cf) == 0) {
          rep(NA, 2)
        }  else if (all(grepl(c("ar1|ar2|ar3|ma1|ma2|ma3|intercept"), names(cf))) &
                    substr(cf["ar1"], 1, j1) %in% arr1 & substr(cf["ar2"], 1, j2) %in% arr2 & substr(cf["ar3"], 1, j3) %in% arr3 & substr(cf["ma1"], 1, k1) %in% maa1 & substr(cf["ma2"], 1, k2) %in% maa2 & substr(cf["ma3"], 1, k3) %in% maa3) {
          c(cf, seed = i)
        } else {
          rep(NA, 2)
        }
      }
      message(' done!\n')

      res1 = res[!sapply(res, anyNA)]


      old <- options()
      on.exit(options(old))
      options(max.print = .Machine$integer.max)

      res2 <- tibble::tibble(Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x)))))

      res2[order(res2$seed), ]

      res2 <- Reduce(function(...) merge(..., all = T), lapply(res1, function(x) as.data.frame(t(x))))
      res2[order(res2$seed), ]
    }

    arima33_sim_search(a = a,  z = z, n = n, p = 3, d = 1, q = 3, ar11 = ar11, ar22 = ar22, ar33 = ar33, ma11 = ma11, ma22 = ma22, ma33 = ma33, sd = sd, j1 = j1, j2 = j2, j3 = j3, k1 = k1, k2 = k2, k3 = k3, arr1 = arr1, arr2 = arr2, arr3 = arr3, maa1 = maa1, maa2 = maa2, maa3 = maa3)

  }
  ##############################################################################

  return(output)
}
