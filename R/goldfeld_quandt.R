#' @description Copied from skedastic::goldfeld_quandt (GPL-3)
#' @source https://github.com/tjfarrar/skedastic/blob/master/R/goldfeld_quandt.R


goldfeld_quandt <- function(mainlm, method = c("parametric", "nonparametric"),
                            deflator = NA, prop_central = 1 / 3, group1prop = 1 / 2,
                            alternative = c("greater", "less", "two.sided"),
                            prob = NA, twosidedmethod = c("doubled", "kulinskaya"),
                            restype = c("ols", "blus"), statonly = FALSE, ...) {
  
  twosidedmethod <- match.arg(twosidedmethod, c("doubled", "kulinskaya"))
  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  method <- match.arg(method, c("parametric", "nonparametric"))
  restype <- match.arg(restype, c("ols", "blus"))
  
  processmainlm(m = mainlm, needy = (method == "parametric"))
  
  hasintercept <- columnof1s(X)
  if (inherits(mainlm, "list")) {
    if (hasintercept[[1]]) {
      if (hasintercept[[2]] != 1) stop("Column of 1's must be first column of design matrix")
      colnames(X) <- c("(Intercept)", paste0("X", 1:(p - 1)))
    } else {
      colnames(X) <- paste0("X", 1:p)
    }
  }
  
  n <- nrow(X)
  
  checkdeflator(deflator, X, p, hasintercept[[1]])
  
  if (method == "parametric") {
    
    theind <- gqind(n, p, prop_central, group1prop)
    
    if (!is.na(deflator) && !is.null(deflator)) {
      if (!is.na(suppressWarnings(as.integer(deflator)))) {
        deflator <- as.integer(deflator)
      }
      y <- y[order(X[, deflator])]
      X <- X[order(X[, deflator]), , drop = FALSE]
    }
    thedf2 <- (length(theind[[2]]) - p)
    thedf1 <- (length(theind[[1]]) - p)
    S2sq <- sum(stats::lm.fit(X[theind[[2]], , drop = FALSE],
                              y[theind[[2]]])$residuals ^ 2) / thedf2
    S1sq <- sum(stats::lm.fit(X[theind[[1]], , drop = FALSE],
                              y[theind[[1]]])$residuals ^ 2) / thedf1
    teststat <- S2sq / S1sq
    if (statonly) return(teststat)
    names(thedf1) <- "df1"
    if (alternative == "greater") {
      pval <- stats::pf(teststat, df1 = thedf1, df2 = thedf2, lower.tail = FALSE)
    } else if (alternative == "less") {
      pval <- stats::pf(teststat, df1 = thedf1, df2 = thedf2, lower.tail = TRUE)
    } else if (alternative == "two.sided") {
      pval <- twosidedpval(q = teststat, CDF = stats::pf,
                           locpar = thedf2 / (thedf2 - 2),
                           method = twosidedmethod, continuous = TRUE,
                           df1 = thedf1, df2 = thedf2,
                           lower.tail = TRUE)
    }
    fullmethod <- "Goldfeld-Quandt F Test"
  } else if (method == "nonparametric") {
    if (restype == "ols") {
      absres <- abs(e)
      newn <- n
    } else if (restype == "blus") {
      absres <- abs(blus(mainlm = list("e" = e, "X" = X), ...))
      newn <- n - p
    }
    
    if (!is.na(deflator) && !is.null(deflator)) {
      if (!is.na(suppressWarnings(as.integer(deflator)))) {
        deflator <- as.integer(deflator)
      }
      absres <- absres[order(X[, deflator])]
    }
    
    teststat <- countpeaks(absres[!is.na(absres)])
    if (statonly) return(teststat)
    
    thedf1 <- NULL
    if (is.na(prob[1]) || is.null(prob)) {
      if (alternative == "greater") {
        pval <- ppeak(k = teststat, n = newn, lower.tail = FALSE, usedata = (newn <= 1000))
      } else if (alternative == "less") {
        pval <- ppeak(k = teststat, n = newn, lower.tail = TRUE, usedata = (newn <= 1000))
      } else if (alternative == "two.sided") {
        peakmean <- sum(0:(newn - 1) * dpeak(k = 0:(newn - 1), n = newn,
                                             usedata = (newn <= 1000)))
        pval <- twosidedpval(q = teststat, locpar = peakmean, CDF = ppeak,
                             method = twosidedmethod, continuous = FALSE,
                             n = newn, lower.tail = TRUE, usedata = (newn <= 1000))
      }
    } else {
      if (length(prob) != newn) stop("prob must be a vector of length equal to number of observations in series")
      if (alternative == "greater") {
        pval <- sum(prob[(teststat + 1):length(prob)])
      } else if (alternative == "less") {
        pval <- sum(prob[1:(teststat + 1)])
      } else if (alternative == "two.sided") {
        pfunc <- function(k) {
          sum(prob[1:(k + 1)])
        }
        peakmean <- sum(0:(newn - 1) * prob)
        pval <- twosidedpval(q = teststat, locpar = peakmean, CDF = pfunc,
                             method = twosidedmethod, continuous = FALSE)
      }
    }
    fullmethod <- "Goldfeld-Quandt Peaks Test"
  } else stop("Invalid `method` argument")
  
  rval <- structure(list(statistic = teststat, p.value = pval, parameter = thedf1,
                         null.value = "Homoskedasticity", alternative = alternative,
                         method = fullmethod), class = "htest")
  broom::tidy(rval)
}