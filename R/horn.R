#' @description Copied from skedastic::horn (GPL-3)
#' @source https://github.com/tjfarrar/skedastic/blob/master/R/horn.R


horn <- function(mainlm, deflator = NA, restype = c("ols", "blus"),
                 alternative = c("two.sided", "greater", "less"),
                 exact = (nres <= 10), statonly = FALSE, ...) {
  
  restype <- match.arg(restype, c("ols", "blus"))
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
  
  processmainlm(m = mainlm, needy = FALSE)
  
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
  
  if (is.na(deflator) || is.null(deflator)) {
    if (restype == "ols") {
      absres <- abs(e)
    } else if (restype == "blus") {
      absres <- abs(blus(mainlm, ...))
      absres <- absres[!is.na(absres)]
    }
  } else {
    if (!is.na(suppressWarnings(as.integer(deflator)))) {
      deflator <- as.integer(deflator)
    }
    if (restype == "ols") {
      absres <- abs(e)[order(X[, deflator])]
    } else if (restype == "blus") {
      absres <- abs(blus(mainlm, ...)[order(X[, deflator])])
      absres <- absres[!is.na(absres)]
    }
  }
  nres <- length(absres)
  
  R <- rank(absres, ties.method = "average")
  teststat <- sum((R - 1:nres) ^ 2)
  if (statonly) return(teststat)
  d <- table(R)
  is.ties <- (max(d) > 1)
  
  if (is.ties) {
    if (exact) {
      warning("Ties are present and exact distribution is not available for D in presence of ties. Normal approximation will be used.")
      exact <- FALSE
    }
    if (max(d) / nres == 1) warning("Normal approximation may not be accurate since maximum rank / m = 1. See Lehmann (1975), p. 294.")
  } else {
    d <- NULL
  }
  
  twosided <- function(tstat, m) {
    if (teststat > (m * (m - 1) * (m + 1) / 6)) {
      return(2 * pDtrend(k = tstat, n = m, lower.tail = FALSE,
                         exact, tiefreq = d))
    } else if (tstat < (m * (m - 1) * (m + 1) / 6)) {
      return(2 * pDtrend(k = tstat, n = m, lower.tail = TRUE,
                         exact, tiefreq = d))
    } else if (tstat == (m * (m - 1) * (m + 1) / 6)) {
      return(1)
    }
  }
  pval <- switch(alternative, "greater" = pDtrend(k = teststat, n = nres,
                                                  lower.tail = FALSE, exact, tiefreq = d),
                 "less" = pDtrend(k = teststat, n = nres,
                                  lower.tail = TRUE, exact, tiefreq = d),
                 "two.sided" = twosided(teststat, m = nres))
  
  rval <- structure(list(statistic = teststat, p.value = pval,
                         alternative = alternative),
                    class = "htest")
  broom::tidy(rval)
}