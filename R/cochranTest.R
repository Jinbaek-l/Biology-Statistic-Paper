#' @description Copied from PMCMRplus::cochranTest (GPL-3)
#' @source https://rdrr.io/cran/PMCMRplus/src/R/cochranTest.R
#' @include pcochran.R


cochranTest <- function(x, ...) UseMethod("cochranTest")

cochranTest.default <- function(x, g, alternative = c("greater", "less"), ...)
{
  if (is.list(x)) {
    if (length(x) < 2L)
      stop("'x' must be a list with at least 2 elements")
    DNAME <- deparse(substitute(x))
    x <- lapply(x, function(u) u <- u[complete.cases(u)])
    k <- length(x)
    l <- sapply(x, "length")
    if (any(l == 0))
      stop("all groups must contain data")
    g <- factor(rep(1 : k, l))
    alternative <- x$alternative
    x <- unlist(x)
  }
  else {
    if (length(x) != length(g))
      stop("'x' and 'g' must have the same length")
    DNAME <- paste(deparse(substitute(x)), "and",
                   deparse(substitute(g)))
    OK <- complete.cases(x, g)
    x <- x[OK]
    g <- g[OK]
    if (!all(is.finite(g)))
      stop("all group levels must be finite")
    g <- factor(g)
    k <- nlevels(g)
    if (k < 2)
      stop("all observations are in the same group")
  }
  
  alternative <- match.arg(alternative)
  Ssq <- tapply(x, g, var)
  n <- length(x)
  
  if (alternative == "greater"){
    method <- "Cochran test for outlying variance"
    TSsq <- max(Ssq)
    i <- which(TSsq == Ssq)
    names(i) <- NULL
    C <- TSsq / sum(Ssq)
    pval <- pcochran(q = C,
                     k = k,
                     n = n / k,
                     lower.tail = FALSE)
    
  } else {
    method <- "Cochran test for inlying variance"
    TSsq <- min(Ssq)
    i <- which(TSsq == Ssq)
    names(i) <- NULL
    C <- TSsq / sum(Ssq)
    pval <- pcochran(q = C,
                     k = k,
                     n = n / k,
                     lower.tail = TRUE)
  }
  
  ans <- list(method = method,
              alternative = alternative,
              statistic = c("C" = C),
              p.value = pval,
              parameter = c("k" = k, "n" = n / k),
              data.name = DNAME,
              estimates = c("group" = i, "var" = TSsq))
  class(ans) <- "htest"
  return(ans)
}

#' @rdname cochranTest
#' @method cochranTest formula
#' @aliases cochranTest.formula
#' @template one-way-formula
#' @export
cochranTest.formula <- function(formula, data, subset, na.action,
                                alternative = c("greater", "less"),...)
{
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  
  if(missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")
  mf <- eval(mf, parent.frame())
  if(length(mf) > 2L)
    stop("'formula' should be of the form response ~ group")
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  alternative <- match.arg(alternative)
  y <- do.call("cochranTest", c(as.list(mf), alternative = alternative))
  y$data.name <- DNAME
  y
}