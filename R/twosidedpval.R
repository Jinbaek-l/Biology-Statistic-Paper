#' @description Copied from skedastic::twosidedpval (GPL-3)
#' @source https://github.com/tjfarrar/skedastic/blob/master/R/twosidedpval.R


twosidedpval <- function(q, CDF, continuous, method = c("doubled", "kulinskaya",
                                                        "minlikelihood"), locpar,
                         supportlim = c(-Inf, Inf), ...) {
  
  method <- match.arg(method, c("doubled", "kulinskaya", "minlikelihood"))
  names(q) <- NULL
  
  if (missing(locpar) && method != "minlikelihood") {
    locpar <- meanfromCDF(theCDF = CDF, cont = continuous, suplim = supportlim,
                          ...)
  }
  
  if (method == "kulinskaya") {
    if (continuous) {
      if (q <= locpar) {
        CDF(q, ...) / CDF(locpar, ...)
      } else {
        (1 - CDF(q, ...)) / (1 - CDF(locpar, ...))
      }
    } else {
      if (!(value_possible(x = locpar, myCDF = CDF, ...))) {
        wL <- CDF(locpar, ...)
        wR <- 1 - CDF(locpar, ...)
        CDF(q, ...) / wL * (q < locpar) + (q == locpar) + (1 - CDF(q - 1, ...)) /
          wR * (q > locpar)
      } else {
        wLm <- CDF(locpar, ...) / (1 + CDF(locpar, ...) - CDF(locpar - 1, ...))
        wRm <- (1 - CDF(locpar - 1, ...)) / (1 + CDF(locpar, ...) - CDF(locpar - 1, ...))
        CDF(q, ...) / wLm * (q < locpar) + (q == locpar) + (1 - CDF(q - 1, ...)) /
          wRm * (q > locpar)
      }
    }
  } else if (method == "doubled") {
    min(2 * CDF(q, ...), 2 * (1 - CDF(q, ...)))
  } else if (method == "minlikelihood") {
    if (continuous) {
      stop("Method minlikelihood not implemented for continuous function")
    } else {
      if (missing(supportlim)) {
        support <- -1e6:1e6
      } else {
        if (supportlim[1] == -Inf) supportlim[1] <- -1e6
        if (supportlim[2] == Inf) supportlim[2] <- 1e6
        support <- supportlim[1]:supportlim[2]
      }
      allPMF <- vapply(support, function(j) CDF(j, ...) -
                         CDF(j - 1, ...), NA_real_)
      PMFq <- CDF(q, ...) - CDF(q - 1, ...)
      sum(allPMF[allPMF <= PMFq])
    }
  } else stop("Invalid `method` argument")
}