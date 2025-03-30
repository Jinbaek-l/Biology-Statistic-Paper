#' @description Copied from skedastic::cook_weisberg (GPL-3)
#' @source https://rdrr.io/cran/skedastic/src/R/cook_weisberg.R


cook_weisberg <- function(mainlm, auxdesign = NA,
                          hetfun = c("mult", "add", "logmult"), statonly = FALSE) {
  
  hetfun <- match.arg(hetfun, c("mult", "add", "logmult"))
  
  auxfitvals <- ifelse(all(is.na(auxdesign)) | is.null(auxdesign), FALSE,
                       auxdesign == "fitted.values")
  processmainlm(m = mainlm, needy = auxfitvals, needyhat = auxfitvals,
                needp = FALSE)
  
  if (all(is.na(auxdesign)) || is.null(auxdesign)) {
    Z <- X
  } else if (is.character(auxdesign)) {
    if (auxdesign == "fitted.values") {
      Z <- t(t(yhat))
    } else stop("Invalid character value for `auxdesign`")
  } else {
    Z <- auxdesign
    if (nrow(auxdesign) != nrow(X)) stop("No. of observations in `auxdesign`
                                         must match\nno. of observations in
                                         original model.")
  }
  
  hasintercept <- columnof1s(Z)
  if (hasintercept[[1]]) {
    Z <- Z[, -hasintercept[[2]], drop = FALSE]
  }
  
  q <- ncol(Z)
  n <- nrow(Z)
  
  if (hetfun == "mult") {
    Z <- cbind(1, Z)
  } else if (hetfun == "logmult") {
    Z <- cbind(1, log(Z))
  } else if (hetfun == "add") {
    Z <- cbind(1, 2 * Z)
  } else stop("Invalid hetfun argument")
  
  sigma_hatsq <- sum(e ^ 2) / n
  std_res_sq <- e ^ 2 / sigma_hatsq
  auxres <- stats::lm.fit(Z, std_res_sq)$residuals
  teststat <- (sum(std_res_sq ^ 2) - n * mean(std_res_sq) ^ 2 - sum(auxres ^ 2)) / 2
  if (statonly) return(teststat)
  method <- hetfun
  pval <- stats::pchisq(teststat, df = q, lower.tail = FALSE)
  rval <- structure(list(statistic = teststat, parameter = q, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "greater", method = method),
                    class = "htest")
  broom::tidy(rval)
}