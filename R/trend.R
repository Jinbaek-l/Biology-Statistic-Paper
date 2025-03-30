#' @description Copied from skedastic::trend (GPL-3)
#' @source https://github.com/tjfarrar/skedastic/blob/master/R/trend.R


dDtrend <- function(k = "all", n, override = FALSE) {
  
  kall <- FALSE
  if (is.na(k[1]) || is.null(k)) {
    stop("Argument k cannot be NA or NULL")
  } else if (is.character(k)) {
    if (length(k) == 1 && k == "all") {
      kall <- TRUE
    } else {
      stop("Invalid character value for k")
    }
  } else if (any(k %% 1 != 0)) {
    stop("Invalid value(s) in k; try an integer or \"all\"")
  }
  
  if (is.na(n[1]) || is.null(n)) {
    stop("Argument n cannot be NA or NULL")
  } else if (!is.numeric(n)) {
    stop("Argument n must be numeric")
  } else if (n %% 1 != 0) {
    stop("Invalid value for n; try an integer")
  } else if (n > 11 && !override) {
    stop("Computation of dDtrend is prohibitively slow for n > 11. Operation aborted. If user insists on proceeding, call function again with `override` set to `TRUE`.")
  }
  
  if (!requireNamespace("arrangements", quietly = TRUE)) {
    stop(
      "Package \"arrangements\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  perms <- arrangements::permutations(n, n)
  valtab <- table(apply(perms, 1, function(x) sum((x - 1:n) ^ 2)))
  prob <- as.double(valtab / factorial(n))
  names(prob) <- names(valtab)
  support <- as.integer(names(prob))
  if (kall) {
    return(prob)
  } else {
    if (any(!(k %in% support))) warning("One or more values in k not part of support of distribution")
    return(prob[support %in% k])
  }
}

pDtrend <- function(k, n, lower.tail = TRUE, exact = (n <= 10), tiefreq = NA,
                    override = FALSE) {
  
  kall <- FALSE
  if (any(is.na(k)) || is.null(k)) {
    stop("Invalid k value(s)")
  } else if (is.character(k)) {
    if (length(k) == 1 && k == "all") {
      kall <- TRUE
    } else {
      stop("Invalid k value(s)")
    }
  } else if (any(!is.numeric(k))) stop("Invalid k value(s)")
  
  if (exact && (is.na(tiefreq[1]) || is.null(tiefreq))) {
    
    prob <- dDtrend(k = "all", n = n, override = override)
    values <- as.integer(names(prob))
    ineqfunc <- ifelse(lower.tail, `<=`, `>=`)
    if (kall) {
      cumprob <- vapply(values, function(j) sum(prob[ineqfunc(values, j)]),
                        NA_real_)
    } else {
      cumprob <- vapply(k, function(j) sum(prob[ineqfunc(values, j)]),
                        NA_real_)
    }
    
  } else {
    
    if (kall) stop("k = \"all\" is only valid when the exact, discrete distribution is used")
    if (is.na(tiefreq[1]) || is.null(tiefreq)) {
      ED <- (n ^ 3 - n) / 6
      VD <- (n ^ 2 * (n + 1) ^ 2 * (n - 1)) / 36
      cumprob <- vapply(k, function(j) stats::pnorm((j - 1 - ED) / sqrt(VD),
                                                    lower.tail = lower.tail), NA_real_)
    } else {
      ED <- (n ^ 3 - n) / 6 - sum(tiefreq ^ 3 - tiefreq) / 12
      VD <- (n ^ 2 * (n + 1) ^ 2 * (n - 1)) / 36 *
        (1 - sum(tiefreq ^ 3 - tiefreq) / (n ^ 3 - n))
      cumprob <- vapply(k, function(j) stats::pnorm((j - ED) / sqrt(VD),
                                                    lower.tail = lower.tail), NA_real_)
    }
  }
  if (kall) {
    names(cumprob) <- names(prob)
  } else {
    names(cumprob) <- as.character(k)
  }
  cumprob
}