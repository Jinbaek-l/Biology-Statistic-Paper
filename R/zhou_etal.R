#' @description Copied from skedastic::zhou_etal (GPL-3)
#' @source https://github.com/tjfarrar/skedastic/blob/master/R/zhou_etal.R


zhou_etal <- function(mainlm, auxdesign = NA,
                      method = c("pooled", "covariate-specific", "hybrid"),
                      Bperturbed = 500L, seed = 1234, statonly = FALSE) {
  
  method <- match.arg(method, c("pooled", "covariate-specific", "hybrid"))
  
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
  if (q == 1 && method %in% c("covariate-specific", "hybrid")) {
    message("Auxiliary design matrix consists of only one column; method changed to `pooled`")
    method <- "pooled"
  }
  n <- nrow(X)
  
  # y <- y - mean(y)
  # X <- apply(X, 2, function(x) x - mean(x))
  # newlm <- stats::lm.fit(X, y)
  # e <- newlm$residuals
  
  if (!is.na(seed)) set.seed(seed)
  Xi <- replicate(Bperturbed, stats::rnorm(n), simplify = FALSE)
  sigmahatsq <- sum(e ^ 2) / n
  H <- Z %*% Rfast::spdinv(crossprod(Z)) %*% t(Z)
  
  if (method == "pooled" || method == "hybrid") {
    wpool <- diag(H) / q
    IRpool <- sum(e ^ 2 * wpool) / sigmahatsq
    Wpoolobs <- sqrt(n) * (IRpool - 1)
    if (statonly && method == "pooled") return(Wpoolobs)
    Wpoolstar <- vapply(Xi, function(xi) sqrt(n) * sum(((wpool - IRpool / n) *
                                                          (e ^ 2 / sigmahatsq - 1) - (IRpool - 1) / n) * xi), NA_real_)
    pvalpool <- sum(Wpoolstar >= Wpoolobs) / Bperturbed
  }
  if (method == "covariate-specific" || method == "hybrid") {
    
    Hminus <- lapply(1:q, function(j) Z[, -j, drop = FALSE] %*%
                       Rfast::spdinv(crossprod(Z[, -j, drop = FALSE])) %*%
                       t(Z[, -j, drop = FALSE]))
    w <- lapply(1:q, function(j) diag(H) - diag(Hminus[[j]]))
    IR <- vapply(1:q, function(j) sum(e ^ 2 * w[[j]]) / sigmahatsq, NA_real_)
    Wobs <- sqrt(n) * (IR - 1)
    if (statonly && method == "covariate-specific") return(Wobs)
    Wstar <- lapply(1:q, function(j) vapply(Xi,
                                            function(xi) sqrt(n) * sum(((w[[j]] - IR[j] / n) *
                                                                          (e ^ 2 / sigmahatsq - 1) - (IR[j] - 1) / n) * xi), NA_real_))
    pvalcov <- q * vapply(1:q, function(j) sum(Wstar[[j]] >= Wobs[j]) /
                            Bperturbed, NA_real_)
    pvalcov[pvalcov > 1] <- 1
    message("covariate-specific p-values include Bonferroni correction")
  }
  
  if (method == "pooled") {
    teststat <- Wpoolobs
    pval <- pvalpool
  } else if (method == "covariate-specific") {
    teststat <- Wobs
    pval <- pvalcov
  } else if (method == "hybrid") {
    minpval <- min(pvalpool, pvalcov)
    pval <- ifelse(2 * minpval > 1, 1, 2 * minpval)
    message("hybrid p-value includes Bonferroni correction")
    if (minpval == pvalpool) {
      teststat <- Wpoolobs
    } else {
      teststat <- Wobs
    }
    if (statonly) return(teststat)
  }
  
  rval <- structure(list(statistic = teststat, parameter = NULL,
                         p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "Heteroskedasticity", method = method), class = "htest")
  broom::tidy(rval)
}