#### Internal utilities: validation, interpolation, quadrature, projection ####
####
#### These helpers turn the (continuous) alpha-level theory of Puri and
#### Ralescu into working, vectorised R code. A fuzzy number is discretised on
#### an alpha grid  0 = a_1 < ... < a_m = 1, with lower envelope
#### a |--> X^-(a) (non-decreasing) and upper envelope a |--> X^+(a)
#### (non-increasing), so that the alpha-cut of X is [X^-(a), X^+(a)].

#' Default alpha grid
#'
#' A fine, regularly spaced grid on the closed unit interval used whenever a
#' level-wise grid is not supplied explicitly by the caller.
#'
#' @return A numeric vector of length 51 on \code{[0, 1]}, containing 0 and 1.
#' @keywords internal
default_alpha <- function() seq(0, 1, length.out = 51L)

#' Trapezoidal quadrature over an alpha grid
#'
#' Approximates \eqn{\int_0^1 f(\alpha)\,d\alpha} by the composite trapezoidal
#' rule applied to the values of \code{f} on the grid \code{alpha}.
#'
#' @param f Numeric vector of function values aligned with \code{alpha}.
#' @param alpha Numeric vector of grid points, sorted in \eqn{[0,1]}.
#' @return A length-one numeric, the approximated integral.
#' @keywords internal
trapz_alpha <- function(f, alpha = default_alpha()) {
  f <- as.numeric(f); alpha <- as.numeric(alpha)
  keep <- is.finite(f) & is.finite(alpha)
  f <- f[keep]; alpha <- alpha[keep]
  if (length(f) < 2L) return(0)
  sum(diff(alpha) * (f[-length(f)] + f[-1L]) / 2)
}

#' Repair a 2x2 covariance matrix for positive (semi)definiteness
#'
#' Symmetrizes \code{Sigma} and, if it is not positive (semi)definite, adds a
#' small diagonal ridge so that a Cholesky factor exists.
#'
#' @param Sigma A 2x2 numeric matrix.
#' @param eps Smallest eigenvalue threshold.
#' @return A symmetric, positive-definite 2x2 numeric matrix.
#' @keywords internal
repair_cov_2x2 <- function(Sigma, eps = 1e-9) {
  Sigma <- 0.5 * (Sigma + t(Sigma))
  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < eps) Sigma <- Sigma + diag(eps * max(1, -min(ev)), 2L)
  Sigma
}

#' Sample from a bivariate normal distribution
#'
#' \code{n} independent draws from \code{N_2(mu, Sigma)} using the Cholesky
#' factorisation (after repair of \code{Sigma}).
#'
#' @param n Number of draws.
#' @param mu Length-2 numeric mean vector.
#' @param Sigma 2x2 covariance matrix.
#' @return An \code{n}x2 numeric matrix; column 1 is the lower endpoint,
#'   column 2 the upper endpoint of a level cut.
#' @keywords internal
mvrnorm_2d <- function(n, mu, Sigma) {
  Sigma <- repair_cov_2x2(Sigma)
  L <- chol(Sigma)
  z <- matrix(stats::rnorm(2L * n), n, 2L)
  sweep(z %*% L, 2L, as.numeric(mu), "+")
}

#' Linearly interpolate the alpha-cuts of a fuzzy number
#'
#' Returns the lower and upper envelopes evaluated on the grid \code{at}.
#'
#' @param alpha Original grid of the fuzzy number (sorted in \code{[0,1]}).
#' @param lower Lower envelope on \code{alpha} (non-decreasing).
#' @param upper Upper envelope on \code{alpha} (non-increasing).
#' @param at Grid on which to evaluate the cuts.
#' @keywords internal
interpolate_cuts <- function(alpha, lower, upper, at) {
  at <- as.numeric(at)
  lf <- stats::approxfun(alpha, lower, rule = 2L,
                         yleft = lower[1L], yright = lower[length(lower)])
  uf <- stats::approxfun(alpha, upper, rule = 2L,
                         yleft = upper[1L], yright = upper[length(upper)])
  list(lower = lf(at), upper = uf(at))
}

#' Interpolate a sample onto a new alpha grid
#'
#' @param obj A \code{FuzzySample} object.
#' @param alpha Target alpha grid.
#' @return An \code{(n, m, 2)} array of cuts on \code{alpha}.
#' @keywords internal
interpolate_sample <- function(obj, alpha) {
  oalpha <- obj$alpha
  n <- obj$n
  m <- length(alpha)
  arr <- array(NA_real_, c(n, m, 2L),
               dimnames = list(NULL, NULL, c("lower", "upper")))
  for (i in seq_len(n)) {
    c <- interpolate_cuts(oalpha, obj$x[i, , 1L], obj$x[i, , 2L], alpha)
    arr[i, , 1L] <- c$lower
    arr[i, , 2L] <- c$upper
  }
  arr
}

#' Project per-level endpoints onto the cone of valid nested cuts
#'
#' Enforces \code{lower <= upper} and monotonicity (lower non-decreasing,
#' upper non-increasing in \code{alpha}).
#'
#' @param lower,upper Numeric vectors, one entry per level.
#' @param alpha Sorted grid in \code{[0,1]}.
#' @return A length-2 list with components \code{lower} and \code{upper}.
#' @keywords internal
project_to_valid_cuts <- function(lower, upper, alpha) {
  lo <- pmin(lower, upper)
  up <- pmax(lower, upper)
  for (k in seq_len(3L)) {
    lo <- stats::isoreg(alpha, lo)$yf
    up <- -stats::isoreg(alpha, -up)$yf
    lo <- pmin(lo, up)
    up <- pmax(lo, up)
  }
  list(lower = lo, upper = up)
}

#' Validate the structure of fuzzy-number cuts
#'
#' @param alpha Sorted grid in \code{[0,1]}.
#' @param lower Non-decreasing lower envelope.
#' @param upper Non-increasing upper envelope, with \code{lower <= upper}.
#' @keywords internal
validate_cuts <- function(alpha, lower, upper) {
  alpha <- as.numeric(alpha); lower <- as.numeric(lower); upper <- as.numeric(upper)
  n <- length(alpha)
  if (n != length(lower) || n != length(upper))
    stop("'alpha', 'lower' and 'upper' must have the same length")
  if (n < 2L) stop("at least two alpha levels are required")
  if (any(diff(alpha) <= 0)) stop("'alpha' must be strictly increasing")
  if (alpha[1L] < 0 || alpha[n] > 1) stop("'alpha' must lie in [0, 1]")
  if (any(lower > upper + 1e-9))
    stop("the lower envelope exceeds the upper envelope at some level")
  if (any(diff(lower) < -1e-9))
    stop("the lower envelope must be non-decreasing in alpha")
    if (any(diff(upper) > 1e-9))
    stop("the upper envelope must be non-increasing in alpha")
  # the core (alpha = 1) is allowed to be an interval [b, c]; only the
  # lower <= upper ordering is enforced (above), so trapezoidal numbers
  # with b < c are valid fuzzy numbers.
  invisible(TRUE)
}

`%||%` <- function(lhs, rhs) if (is.null(lhs)) rhs else lhs

#' Standardise an alpha grid
#'
#' Ensures a strictly increasing grid on \code{[0, 1]} containing 0 and 1.
#'
#' @param alpha User-supplied alpha grid (may be \code{NULL}).
#' @keywords internal
standardise_alpha <- function(alpha = NULL) {
  if (is.null(alpha)) return(default_alpha())
  alpha <- as.numeric(alpha)
  if (anyDuplicated(alpha)) stop("duplicate alpha levels are not allowed")
  alpha <- sort(alpha)
  if (alpha[1L] < 0 || alpha[length(alpha)] > 1)
    stop("alpha levels must lie in [0, 1]")
  if (alpha[1L] > 0) alpha <- c(0, alpha)
  if (alpha[length(alpha)] < 1) alpha <- c(alpha, 1)
  alpha
}
