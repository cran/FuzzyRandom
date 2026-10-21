#### Fuzzy numbers (alpha-cut representation) ####
####
#### A fuzzy number A is stored as the triple (alpha, lower, upper) where
#### alpha in [0, 1] and the alpha-cut of A is [A^-(alpha), A^+(alpha)] =
#### [lower, upper], with lower non-decreasing and upper non-increasing in
#### alpha. This is the representation used throughout Puri-Ralescu theory.

#' A fuzzy number given by its alpha-level cuts
#'
#' Build a \code{\link[=FuzzyNumber]{FuzzyNumber}} object from its
#' alpha-level envelopes. When \code{alpha} is not supplied it is taken as a
#' regularly spaced grid on \code{[0, 1]} with one node per supplied endpoint.
#'
#' @param alpha Numeric vector in \code{[0, 1]}, strictly increasing.
#' @param lower Numeric; lower envelope, non-decreasing in \code{alpha}.
#' @param upper Numeric; upper envelope, non-increasing in \code{alpha}.
#' @return An object of class \code{"FuzzyNumber"}.
#' @examples
#' \dontshow{knitr::opts_chunk$set(exdent = TRUE)}
#' fn <- FuzzyNumber(alpha = c(0, 0.5, 1), lower = c(1, 2, 2), upper = c(4, 2, 2))
#' cuts(fn)
#' @export
FuzzyNumber <- function(alpha = NULL, lower, upper) {
  if (missing(lower) || missing(upper))
    stop("please supply 'lower' and 'upper'")
  if (length(lower) != length(upper))
    stop("'lower' and 'upper' must have the same length")
  if (is.null(alpha))
    alpha <- seq(0, 1, length.out = length(lower))
  alpha <- standardise_alpha(alpha)
  if (length(lower) != length(alpha))
    stop("lengths of 'alpha', 'lower' and 'upper' must match")
  validate_cuts(alpha, lower, upper)
  structure(list(alpha = as.numeric(alpha),
                 lower = as.numeric(lower),
                 upper = as.numeric(upper)),
            class = "FuzzyNumber")
}

#' Triangular fuzzy number
#'
#' The triangular fuzzy number with support \code{[a, c]} and core \code{b}
#' (\code{a <= b <= c}). Its alpha-cuts are
#' \code{[a + (b-a)*alpha, c - (c-b)*alpha]}.
#'
#' @param a,b,c Scalars with \code{a <= b <= c}.
#' @param alpha Optional alpha grid (see \code{\link{FuzzyNumber}}).
#' @return An object of class \code{"FuzzyNumber"}.
#' @examples
#' plot(triangular_fuzzy(1, 2, 4))
#' @export
triangular_fuzzy <- function(a, b, c, alpha = NULL) {
  a <- as.numeric(a); b <- as.numeric(b); c <- as.numeric(c)
  if (length(a) != 1L || length(b) != 1L || length(c) != 1L)
    stop("'a', 'b' and 'c' must be scalars")
  if (!(a <= b && b <= c)) stop("triangular_fuzzy requires a <= b <= c")
  alpha <- standardise_alpha(alpha)
  FuzzyNumber(alpha, a + (b - a) * alpha, c - (c - b) * alpha)
}

#' Trapezoidal fuzzy number
#'
#' The trapezoidal fuzzy number with support \code{[a, d]} and core
#' \code{[b, c]} (\code{a <= b <= c <= d}).
#'
#' @param a,b,c,d Scalars with \code{a <= b <= c <= d}.
#' @param alpha Optional alpha grid.
#' @return An object of class \code{"FuzzyNumber"}.
#' @examples
#' plot(trapezoidal_fuzzy(1, 2, 3, 5))
#' @export
trapezoidal_fuzzy <- function(a, b, c, d, alpha = NULL) {
  a <- as.numeric(a); b <- as.numeric(b); c <- as.numeric(c); d <- as.numeric(d)
  if (length(a) != 1L || length(b) != 1L || length(c) != 1L || length(d) != 1L)
    stop("'a', 'b', 'c' and 'd' must be scalars")
  if (!(a <= b && b <= c && c <= d))
    stop("trapezoidal_fuzzy requires a <= b <= c <= d")
  alpha <- standardise_alpha(alpha)
  FuzzyNumber(alpha, a + (b - a) * alpha, d - (d - c) * alpha)
}

#' A fuzzy number from a membership function evaluated on a grid
#'
#' Reconstructs the alpha-cuts of a (normal, quasi-concave) fuzzy number from
#' a membership function evaluated on a grid. The level set
#' \code{\{x : mu(x) >= alpha\}} is assumed to be a single interval.
#'
#' @param x Sorted numeric grid of the support (unsorted grids are sorted).
#' @param mu Membership values in \code{[0, 1]} evaluated at \code{x}.
#' @param alpha Optional alpha grid.
#' @return An object of class \code{"FuzzyNumber"}.
#' @examples
#' x  <- seq(-3, 3, length.out = 201)
#' mu <- pmax(1 - abs(x), 0)            # triangular membership
#' fn <- fuzzy_from_membership(x, mu)
#' plot(fn)
#' @export
fuzzy_from_membership <- function(x, mu, alpha = NULL) {
  if (length(x) != length(mu))
    stop("'x' and 'mu' must have the same length")
  if (!is.numeric(x) || !is.numeric(mu))
    stop("'x' and 'mu' must be numeric")
  if (!is.unsorted(x)) { } else { o <- sort.list(x); x <- x[o]; mu <- mu[o] }
  mu <- pmax(pmin(mu, 1), 0)
  alpha <- standardise_alpha(alpha)
  peak <- which.max(mu)
  lx <- x[seq_len(peak)];     lmu <- mu[seq_len(peak)]      # increasing
  rx <- x[peak:length(x)];    rmu <- mu[peak:length(x)]     # decreasing
  lower <- stats::approxfun(lmu, lx, rule = 2L,
                            yleft = lx[1L], yright = lx[peak])(alpha)
  upper <- stats::approxfun(rev(rmu), rev(rx), rule = 2L,
                            yleft = rx[length(rx)], yright = rx[1L])(alpha)
  lower <- pmin(lower, upper)
  upper <- pmax(lower, upper)
  FuzzyNumber(alpha, lower, upper)
}

#### Cuts, accessors and methods for FuzzyNumber ####

#' Alpha-level cuts of a fuzzy number
#'
#' @param x A `FuzzyNumber` (or other compatible object).
#' @param alpha Optional alpha grid at which to evaluate the cuts.
#' @param ... Unused.
#' @return A numeric matrix with columns `lower` and `upper`.
#' @examples
#' fn <- triangular_fuzzy(1, 2, 4)
#' cuts(fn, alpha = c(0, 0.5, 1))
#' @rdname cuts
#' @export
cuts <- function(x, alpha = NULL, ...) {
  if (is.null(alpha)) alpha <- default_alpha()
  UseMethod("cuts")
}

#' @method cuts FuzzyNumber
#' @rdname cuts
#' @export
cuts.FuzzyNumber <- function(x, alpha = NULL, ...) {
  if (is.null(alpha)) alpha <- x$alpha
  c <- interpolate_cuts(x$alpha, x$lower, x$upper, as.numeric(alpha))
  matrix(c(c$lower, c$upper), length(alpha), 2L,
         dimnames = list(format(alpha), c("lower", "upper")))
}

#' Cut of a fuzzy number at a single level
#'
#' @param x A `FuzzyNumber`.
#' @param alpha A single level in `[0, 1]`.
#' @param ... Unused.
#' @return A named numeric vector `c(lower = ..., upper = ...)`.
#' @examples
#' cut_at(triangular_fuzzy(1, 2, 4), 0.5)
#' @export
cut_at <- function(x, alpha, ...) UseMethod("cut_at")

#' @method cut_at FuzzyNumber
#' @rdname cut_at
#' @export
cut_at.FuzzyNumber <- function(x, alpha = 0.5, ...) {
  if (length(alpha) != 1L) stop("'alpha' must be a single level")
  c <- interpolate_cuts(x$alpha, x$lower, x$upper, as.numeric(alpha))
  c(lower = c$lower[1L], upper = c$upper[1L])
}

#' Number of alpha levels of a fuzzy object
#'
#' @param x A `FuzzyNumber` or `FuzzySample`.
#' @return Length-one integer.
#' @examples
#' n_cuts(triangular_fuzzy(1, 2, 4))
#' @export
n_cuts <- function(x) {
  if (inherits(x, "FuzzyNumber")) return(length(x$alpha))
  if (inherits(x, "FuzzySample")) return(length(x$alpha))
  stop("n_cuts requires a FuzzyNumber or a FuzzySample")
}

#' Alpha grid of a fuzzy object
#'
#' @param x A `FuzzyNumber` or `FuzzySample`.
#' @return Numeric vector of alpha levels.
#' @examples
#' alpha_levels(triangular_fuzzy(1, 2, 4))[1:3]
#' @export
alpha_levels <- function(x) {
  if (inherits(x, "FuzzyNumber")) return(x$alpha)
  if (inherits(x, "FuzzySample")) return(x$alpha)
  stop("alpha_levels requires a FuzzyNumber or a FuzzySample")
}

#' Whether an object is a (valid) fuzzy number
#'
#' @param x Any R object.
#' @return Logical; `TRUE` if `x` is a valid `FuzzyNumber`.
#' @examples
#' is_fuzzy_number(triangular_fuzzy(0, 1, 2))
#' is_fuzzy_number(1:3)
#' @export
is_fuzzy_number <- function(x) {
  inherits(x, "FuzzyNumber") &&
    is.numeric(x$alpha) && is.numeric(x$lower) && is.numeric(x$upper) &&
    tryCatch(validate_cuts(x$alpha, x$lower, x$upper),
             error = function(e) FALSE)
}

#' @method print FuzzyNumber
#' @export
print.FuzzyNumber <- function(x, ...) {
  cat("<FuzzyNumber>\n")
  cat(sprintf("  levels: %d  (alpha %g .. %g)\n",
              length(x$alpha), x$alpha[1L], x$alpha[length(x$alpha)]))
  cat(sprintf("  support : [ %g , %g ]\n", x$lower[1L], x$upper[1L]))
  n <- length(x$alpha)
  cat(sprintf("  core    : [ %g , %g ]\n", x$lower[n], x$upper[n]))
  invisible(x)
}

#' @method summary FuzzyNumber
#' @export
summary.FuzzyNumber <- function(object, ...) {
  alpha <- object$alpha
  n <- length(alpha)
  structure(
    list(n_levels = n,
         support = c(object$lower[1L], object$upper[1L]),
         core = c(object$lower[n], object$upper[n]),
         mid_spread = cbind((object$lower + object$upper) / 2,
                            (object$upper - object$lower) / 2)),
    class = "summary_FuzzyNumber")
}

#' @method print summary_FuzzyNumber
#' @export
print.summary_FuzzyNumber <- function(x, ...) {
  cat("<summary FuzzyNumber>\n")
  cat(sprintf("  levels : %d\n  support: [ %g , %g ]\n  core   : [ %g , %g ]\n",
              x$n_levels, x$support[1L], x$support[2L], x$core[1L], x$core[2L]))
  invisible(x)
}

#' @method as.data.frame FuzzyNumber
#' @export
as.data.frame.FuzzyNumber <- function(x, ...) {
  data.frame(alpha = x$alpha, lower = x$lower, upper = x$upper)
}

#' @method plot FuzzyNumber
#' @export
plot.FuzzyNumber <- function(x, main = "Fuzzy number",
                             value_lab = "value", col = "steelblue",
                             border = NA, ...) {
  alpha <- x$alpha; lo <- x$lower; up <- x$upper
  xlim <- range(c(lo, up))
  graphics::plot(xlim, c(0, 1), type = "n", xlab = value_lab,
                 ylab = "alpha", main = main, xlim = xlim, ylim = c(0, 1), ...)
  graphics::polygon(c(lo, rev(up)), c(alpha, rev(alpha)),
                    col = grDevices::rgb(0.26, 0.49, 0.85, 0.25), border = border)
  graphics::lines(lo, alpha, col = col, lwd = 2)
  graphics::lines(up, alpha, col = col, lwd = 2)
  invisible(x)
}

