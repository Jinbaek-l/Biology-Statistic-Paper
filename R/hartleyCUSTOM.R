#' @description Adapted from PMCMRplus::hartleyTest (GPL-3)
#' @source https://rdrr.io/cran/PMCMRplus/src/R/hartleyTest.R
#' @include SuppDists package (call an externally function implemented by C) 

#' @details Error when the minimum variance, or denominator becomes zero. For this case, set the minimum variance to 0.0002675227, which was empirically obtained from COPDGene project data.

hartleyCUSTOM <- function(x, ...) UseMethod("hartleyCUSTOM")

#' @rdname hartleyCUSTOM
#' @method hartleyCUSTOM default
#' @template one-way-parms
#' @importFrom SuppDists pmaxFratio
#' @export
hartleyCUSTOM.default <-
  function(x,
           g,
           ...)
  {
    ## taken from stats::hartleyTest
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
    
    
    ##
    Var <- tapply(x, g, var)
    n <- tapply(x, g, length)
    N <- mean(n)
    k <- nlevels(g)
    
    if (any(n != N)) {
      warning("Maximum F-ratio test is imprecise for unbalanced designs.")
    }
    mxId <- which.max(Var)
    mnId <- which.min(Var)
    
    ## Choose n of minimum variance,
    ## if n is equal, than thiswill not have any effect
    df <- n[mnId] - 1
    
    if(Var[mnId]==0){
      ### Code modified!!!!! 
      PSTAT <- Var[mxId] / 0.0002675227 
    }else{
      PSTAT <- Var[mxId] / Var[mnId]
    }
    
    names(PSTAT) <- "F Max"
    PARMS <- c(df, k)
    names(PARMS) <- c("df", "k")
    
    PVAL <- pmaxFratio(PSTAT, df = df, k = k, lower.tail = FALSE)
    
    METHOD <-
      paste("Hartley's maximum F-ratio test of homogeneity of variances")
    
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, parameter = PARMS)
    class(ans) <- "htest"
    ans
  }

#' @rdname hartleyCUSTOM
#' @method hartleyCUSTOM formula
#' @template one-way-formula
#' @export
hartleyCUSTOM.formula <-
  function(formula, data, subset, na.action,
           ...)
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
    y <- do.call("hartleyCUSTOM", c(as.list(mf)))
    y$data.name <- DNAME
    y
  }