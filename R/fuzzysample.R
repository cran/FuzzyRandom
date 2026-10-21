#### Samples of fuzzy random variables (FuzzySample) ####
####
#### A sample X_1, ..., X_n of i.i.d. fuzzy random variables is stored as an
#### (n, m, 2) array `x` of alpha-cuts, the third dimension being
#### c("lower", "upper"), together with the shared alpha grid.

#' Build a sample of fuzzy random variables
#'
#' @param x Either a list of \code{FuzzyNumber} objects, an \code{(n, m, 2)}
#'   numeric array of cuts (dimension 3 = \code{c("lower", "upper")}), or an
#'   existing \code{FuzzySample}.
#' @param alpha Alpha grid shared by all observations. Required when \code{x}
#'   is an array with an unlabelled third dimension; otherwise inferred.
#' @return An object of class \code{"FuzzySample"}.
#' @examples
#' set.seed(1)
#' s <- FuzzySample(list(
#'   triangular_fuzzy(0, 1, 2),
#'   triangular_fuzzy(1, 2, 3),
#'   trapezoidal_fuzzy(0.5, 1, 3, 4)))
#' n_obs(s)
#' @export
FuzzySample <- function(x, alpha = NULL) {
  if (inherits(x, "FuzzySample")) {
    if (!is.null(alpha)) {
      a <- standardise_alpha(alpha)
      arr <- interpolate_sample(x, a)
      return(new_FuzzySample(arr, a))
    }
    return(x)
  }
  if (is.list(x) && length(x) > 0L && inherits(x[[1L]], "FuzzyNumber")) {
    a <- standardise_alpha(alpha)
    n <- length(x); m <- length(a)
    arr <- array(NA_real_, c(n, m, 2L),
                 dimnames = list(NULL, NULL, c("lower", "upper")))
    for (i in seq_len(n)) {
      c <- interpolate_cuts(x[[i]]$alpha, x[[i]]$lower,
                            x[[i]]$upper, a)
      arr[i, , 1L] <- c$lower
      arr[i, , 2L] <- c$upper
    }
    for (i in seq_len(n)) validate_cuts(a, arr[i, , 1L], arr[i, , 2L])
    return(new_FuzzySample(arr, a))
  }
  if (is.array(x) && length(dim(x)) == 3L) {
    d <- dim(x)
    if (is.null(alpha)) alpha <- seq(0, 1, length.out = d[2L])
    a <- standardise_alpha(alpha)
    if (length(a) != d[2L])
      stop("length of 'alpha' must match dim 2 of the array")
    if (is.null(dimnames(x)[[3L]])) dimnames(x)[[3L]] <- c("lower", "upper")
    return(new_FuzzySample(x, a))
  }
  stop("cannot coerce 'x' to a FuzzySample; supply a list of FuzzyNumber or an ",
       "(n, m, 2) array of cuts")
}

#' Internal constructor (no validation) for FuzzySample objects
#'
#' @param arr An \code{(n, m, 2)} array of cuts.
#' @param alpha The shared alpha grid.
#' @return A \code{FuzzySample} object.
#' @keywords internal
new_FuzzySample <- function(arr, alpha) {
  structure(list(x = arr, alpha = as.numeric(alpha), n = dim(arr)[1L]),
            class = "FuzzySample")
}

#' Coerce to a FuzzySample
#'
#' Coerces from a list of \code{FuzzyNumber}, an array, a data frame (with
#' columns `id`, `alpha`, `lower`, `upper`) or a single \code{FuzzyNumber}
#' (giving a sample of size 1).
#'
#' @param x Object to coerce.
#' @param alpha Alpha grid (for array/data-frame input).
#' @param ... Unused.
#' @return A \code{FuzzySample}.
#' @examples
#' fs <- as_fuzzy_sample(triangular_fuzzy(0, 1, 2))
#' n_obs(fs)
#' @rdname as_fuzzy_sample
#' @export
as_fuzzy_sample <- function(x, alpha = NULL, ...) {
  if (inherits(x, "FuzzySample")) {
    if (!is.null(alpha)) return(FuzzySample(x, alpha))
    return(x)
  }
  UseMethod("as_fuzzy_sample")
}

#' @method as_fuzzy_sample default
#' @rdname as_fuzzy_sample
#' @export
as_fuzzy_sample.default <- function(x, alpha = NULL, ...) {
  FuzzySample(x, alpha)
}

#' @method as_fuzzy_sample FuzzyNumber
#' @rdname as_fuzzy_sample
#' @export
as_fuzzy_sample.FuzzyNumber <- function(x, alpha = NULL, ...) {
  FuzzySample(list(x), alpha)
}

#' @method as_fuzzy_sample data.frame
#' @rdname as_fuzzy_sample
#' @export
as_fuzzy_sample.data.frame <- function(x, alpha = NULL, ...) {
  req <- c("alpha", "lower", "upper")
  if (!all(req %in% names(x)))
    stop("data frame must contain 'alpha', 'lower' and 'upper'")
  if (!is.null(x$id)) {
    ids <- sort(unique(x$id))
    a <- standardise_alpha(alpha)
    n <- length(ids); m <- length(a)
    arr <- array(NA_real_, c(n, m, 2L),
                 dimnames = list(NULL, NULL, c("lower", "upper")))
    for (k in seq_len(n)) {
      rowk <- x$id == ids[k]
      xk <- x[rowk, , drop = FALSE]
      ak <- sort(xk$alpha)
      arr[k, , 1L] <- stats::approx(ak, xk$lower[order(xk$alpha)],
                                    xout = a, rule = 2L)$y
      arr[k, , 2L] <- stats::approx(ak, xk$upper[order(xk$alpha)],
                                    xout = a, rule = 2L)$y
    }
    return(new_FuzzySample(arr, a))
  }
  FuzzySample(as.matrix(x[, c("lower", "upper")]), x$alpha)
}


#### FuzzySample accessors and methods ####

#' Number of observations in a sample
#'
#' @param x A \code{FuzzySample}.
#' @return Length-one integer.
#' @examples
#' s <- FuzzySample(list(triangular_fuzzy(0, 1, 2), triangular_fuzzy(1, 2, 3)))
#' n_obs(s)
#' @export
n_obs <- function(x) {
  if (inherits(x, "FuzzySample")) return(x$n)
  stop("n_obs requires a FuzzySample")
}

#' Whether an object is a FuzzySample
#'
#' @param x Any R object.
#' @return Logical.
#' @examples
#' is_fuzzy_sample(FuzzySample(list(triangular_fuzzy(0, 1, 2))))
#' @export
is_fuzzy_sample <- function(x) inherits(x, "FuzzySample")

#' @method print FuzzySample
#' @export
print.FuzzySample <- function(x, ...) {
  cat("<FuzzySample>\n")
  cat(sprintf("  n    : %d fuzzy random variables\n", x$n))
  cat(sprintf("  cuts : %d alpha levels (alpha %g .. %g)\n",
              length(x$alpha), x$alpha[1L], x$alpha[length(x$alpha)]))
  invisible(x)
}

#' @method summary FuzzySample
#' @export
summary.FuzzySample <- function(object, ...) {
  a <- object$alpha
  n <- object$n
  lo <- object$x[, , 1L]; up <- object$x[, , 2L]
  structure(
    list(n = n, n_levels = length(a),
         support = matrix(range(lo), 2L, 2L,
                          dimnames = list(c("min", "max"), c("lower", "upper"))),
         mean = cbind(colMeans(lo), colMeans(up))),
    class = "summary_FuzzySample")
}

#' @method print summary_FuzzySample
#' @export
print.summary_FuzzySample <- function(x, ...) {
  cat("<summary FuzzySample>\n")
  cat(sprintf("  n          : %d\n", x$n))
  cat(sprintf("  levels     : %d\n", x$n_levels))
  cat(sprintf("  support    : lower [%g , %g], upper [%g , %g]\n",
              x$support["min", "lower"], x$support["max", "lower"],
              x$support["min", "upper"], x$support["max", "upper"]))
  cat("  mean cut   : lower / upper (per level):\n")
  print.default(x$mean)
  invisible(x)
}

#' @method as.data.frame FuzzySample
#' @export
as.data.frame.FuzzySample <- function(x, ..., id = "id") {
  a <- x$alpha; n <- x$n; arr <- x$x
  df <- data.frame(
    id = rep(seq_len(n), each = length(a)),
    alpha = rep(a, n),
    lower = as.vector(arr[, , 1L]),
    upper = as.vector(arr[, , 2L]))
  names(df)[1L] <- id
  df
}

#' @method print as.data.frame.FuzzySample
#' @export
print.as.data.frame.FuzzySample <- function(x, ...) {
  print.data.frame(x, ...)
}

#' @method [ FuzzySample
#' @export
'[.FuzzySample' <- function(x, i, ...) {
  arr <- x$x[i, , , drop = FALSE]
  if (is.null(dim(arr))) {
    arr <- array(arr, c(length(i), length(x$alpha), 2L),
                 dimnames = list(NULL, NULL, c("lower", "upper")))
  }
  new_FuzzySample(arr, x$alpha)
}

#' @method plot FuzzySample
#' @export
plot.FuzzySample <- function(x, main = "Fuzzy sample",
                             value_lab = "value",
                             max_obs = 50, show_mean = TRUE, ...) {
  a <- x$alpha; arr <- x$x
  obs <- if (x$n <= max_obs) seq_len(x$n) else sample(x$n, max_obs)
  xlim <- range(c(arr[obs, , 1L], arr[obs, , 2L]))
  graphics::plot(xlim, c(0, 1), type = "n", xlab = value_lab,
                 ylab = "alpha", main = main, xlim = xlim, ylim = c(0, 1), ...)
  for (k in obs) {
    lo <- arr[k, , 1L]; up <- arr[k, , 2L]
    graphics::polygon(c(lo, rev(up)), c(a, rev(a)),
                      col = grDevices::rgb(0.26, 0.49, 0.85, 0.05), border = NA)
  }
  if (show_mean) {
    ml <- colMeans(arr[, , 1L]); mu <- colMeans(arr[, , 2L])
    graphics::polygon(c(ml, rev(mu)), c(a, rev(a)),
                      col = grDevices::rgb(0.86, 0.15, 0.16, 0.35), border = NA)
    graphics::lines(ml, a, col = "red", lwd = 2)
    graphics::lines(mu, a, col = "red", lwd = 2)
  }
  invisible(x)
}
