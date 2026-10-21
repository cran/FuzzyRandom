#### Expectation, variance, covariance (Puri-Ralescu, level-wise) ####
####
#### The expectation is E[X]_alpha = [E[X^-(alpha)], E[X^+(alpha)]];
#### covariance of X,Y at level alpha is the 2x2 matrix
####   [[Cov(X^-,Y^-), Cov(X^-,Y^+)],
####    [Cov(X^+,Y^-), Cov(X^+,Y^+)]];
#### the variance is Cov(X,X); the scalar reductions integrate the lower and
#### upper diagonal entries over alpha (Feng--Hu--Shu, 2002).

#' Internal constructor for a FuzzyNumber (no validation)
#'
#' @param alpha Alpha grid.
#' @param lower,upper Envelopes.
#' @keywords internal
new_FuzzyNumber <- function(alpha, lower, upper)
  structure(list(alpha = as.numeric(alpha),
                 lower = as.numeric(lower),
                 upper = as.numeric(upper)),
            class = "FuzzyNumber")

#' Internal constructor for a FuzzyCovariance object
#'
#' @param alpha Alpha grid.
#' @param cov Array of dim \code{c(m, 2, 2)} (rows/cols: lower, upper).
#' @param scalar Scalar reduction (Feng--Hu--Shu integral).
#' @param type \code{"variance"} or \code{"covariance"}.
#' @keywords internal
new_FuzzyCovariance <- function(alpha, cov, scalar, type) {
  structure(
    list(alpha = as.numeric(alpha), cov = cov,
         var_lower = cov[, 1L, 1L], var_upper = cov[, 2L, 2L],
         cov_cross = cov[, 1L, 2L],
         scalar = scalar, type = type),
    class = "FuzzyCovariance")
}

#' Expectation of a fuzzy random variable (sample)
#'
#' Level-wise expectation.  For a sample the estimator is the level-wise mean
#' of the alpha-cuts:
#' \deqn{[E[X]]_{\alpha} = \Bigl[\,\frac1n\sum_i X_i^{-}(\alpha),\;
#'       \frac1n\sum_i X_i^{+}(\alpha)\Bigr].}
#'
#' @param x A \code{FuzzySample}.
#' @param alpha Optional alpha grid on which to evaluate.
#' @param ... Unused.
#' @return A \code{FuzzyNumber}, the sample expectation.
#' @examples
#' set.seed(2)
#' s <- FuzzySample(replicate(40, triangular_fuzzy(-1, 0, 1), simplify = FALSE))
#' ex <- expectation(s); plot(ex, main = "sample mean")
#' @export
expectation <- function(x, alpha = NULL, ...) {
  if (is.null(alpha)) alpha <- x$alpha
  UseMethod("expectation")
}

#' @method expectation FuzzySample
#' @rdname expectation
#' @export
expectation.FuzzySample <- function(x, alpha = NULL, ...) {
  if (is.null(alpha)) alpha <- x$alpha
  alpha <- standardise_alpha(alpha)
  arr <- interpolate_sample(x, alpha)
  lo <- colMeans(arr[, , 1L])
  up <- colMeans(arr[, , 2L])
  new_FuzzyNumber(alpha, lo, up)
}

#' @method mean FuzzySample
#' @rdname expectation
#' @export
mean.FuzzySample <- function(x, ...) expectation(x, ...)

#' Covariance of two fuzzy random variables (sample)
#'
#' Level-wise covariance.  At each level \code{alpha} the covariance of
#' \code{X} and \code{Y} is the \code{2 x 2} matrix
#' \tabular{ll}{
#' \code{Cov(X^-, Y^-)} \tab \code{Cov(X^-, Y^+)} \cr
#' \code{Cov(X^+, Y^-)} \tab \code{Cov(X^+, Y^+)}.
#' }
#' When \code{y} is omitted, \code{x} is paired with itself and the result is
#' the variance.  A scalar reduction through the Feng--Hu--Shu integral
#' \deqn{\mathrm{Cov}_{FHS}(X,Y)=\tfrac12\int_0^1\bigl[\mathrm{Cov}(X^-(a),Y^-(a))+
#'        \mathrm{Cov}(X^+(a),Y^+(a))\bigr]\,da} is also returned.
#'
#' @param x,y Two \code{FuzzySample} objects (or \code{x} alone for the
#'   variance).
#' @param alpha Optional alpha grid.
#' @param ... Unused.
#' @return A \code{FuzzyCovariance} object.
#' @examples
#' set.seed(3)
#' x <- FuzzySample(replicate(30, triangular_fuzzy(0, 1, 2), simplify = FALSE))
#' y <- FuzzySample(replicate(30, triangular_fuzzy(1, 2, 3), simplify = FALSE))
#' covariance(x, y)
#' @export
covariance <- function(x, y = NULL, alpha = NULL, ...) {
  if (is.null(y)) alpha <- x$alpha
  UseMethod("covariance")
}

#' @method covariance FuzzySample
#' @rdname covariance
#' @export
covariance.FuzzySample <- function(x, y = NULL, alpha = NULL, ...) {
  if (is.null(y)) y <- x
  if (is.null(alpha)) alpha <- x$alpha
  alpha <- standardise_alpha(alpha)
  xa <- interpolate_sample(x, alpha)
  ya <- interpolate_sample(y, alpha)
  n <- x$n; m <- length(alpha)
  cov_arr <- array(NA_real_, c(m, 2L, 2L),
                   dimnames = list(NULL, c("lower", "upper"), c("lower", "upper")))
  if (n > 1L) {
    for (j in seq_len(m)) {
      A <- matrix(xa[, j, ], n, 2L)
      B <- matrix(ya[, j, ], n, 2L)
      cov_arr[j, , ] <- stats::cov(A, B)
    }
  } else {
    warning("covariance is not defined for a sample of size 1")
  }
  scalar <- 0.5 * trapz_alpha(cov_arr[, 1L, 1L] + cov_arr[, 2L, 2L], alpha)
  type <- if (isTRUE(identical(y, x))) "variance" else "covariance"
  new_FuzzyCovariance(alpha, cov_arr, scalar, type)
}

#' @method cov FuzzySample
#' @rdname covariance
#' @export
cov.FuzzySample <- function(x, y = NULL, ...) covariance(x, y, ...)

#' Variance of a fuzzy random variable (sample)
#'
#' Level-wise variance, i.e. \code{\link{covariance}} of a sample with itself.
#' Returns the level-wise \code{2 x 2} matrix
#' \code{[[Var(X^-), Cov(X^-,X^+)], [Cov(X^+,X^-), Var(X^+)]} at each level,
#' together with the scalar Feng--Hu--Shu variance
#' \deqn{\mathrm{Var}_{FHS}(X)=\tfrac12\int_0^1\bigl[\mathrm{Var}(X^-(a))+
#'        \mathrm{Var}(X^+(a))\bigr]\,da.}
#'
#' @param x A \code{FuzzySample}.
#' @param alpha Optional alpha grid.
#' @param ... Unused.
#' @return A \code{FuzzyCovariance} object (with \code{type = "variance"}).
#' @examples
#' set.seed(4)
#' s <- FuzzySample(replicate(30, triangular_fuzzy(0, 1, 2), simplify = FALSE))
#' variance(s)
#' @export
variance <- function(x, alpha = NULL, ...) {
  UseMethod("variance")
}

#' @method variance FuzzySample
#' @rdname variance
#' @export
variance.FuzzySample <- function(x, alpha = NULL, ...) covariance(x, x, alpha)

#' @method var FuzzySample
#' @rdname variance
#' @export
var.FuzzySample <- function(x, ...) variance(x, ...)

#### Correlation, scalar (Feng--Hu--Shu) reductions, and printing ####

#' Correlation of two fuzzy random variables (sample)
#'
#' Level-wise Pearson correlation between the lower endpoints
#' (\code{cor(X^-(alpha), Y^-(alpha))}) and between the upper endpoints
#' (\code{cor(X^+(alpha), Y^+(alpha))}), plus the scalar reduction
#' \deqn{\rho_{FHS}(X,Y)=\tfrac12\int_0^1\bigl[r_{--}(a)+r_{++}(a)\bigr]\,da.}
#'
#' @param x,y Two \code{FuzzySample} objects.
#' @param alpha Optional alpha grid.
#' @param ... Unused.
#' @return An object of class \code{"FuzzyCorrelation"} (with a scalar slot).
#' @examples
#' set.seed(5)
#' x <- FuzzySample(replicate(50, triangular_fuzzy(0, 1, 2), simplify = FALSE))
#' y <- FuzzySample(replicate(50, triangular_fuzzy(0, 1, 2), simplify = FALSE))
#' correlation(x, y)
#' @export
correlation <- function(x, y = NULL, alpha = NULL, ...) {
  if (is.null(y)) y <- x
  UseMethod("correlation")
}

#' @method correlation FuzzySample
#' @rdname correlation
#' @export
correlation.FuzzySample <- function(x, y = NULL, alpha = NULL, ...) {
  if (is.null(y)) y <- x
  if (is.null(alpha)) alpha <- x$alpha
  alpha <- standardise_alpha(alpha)
  xa <- interpolate_sample(x, alpha)
  ya <- interpolate_sample(y, alpha)
  n <- x$n; m <- length(alpha)
  rll <- rep(NA_real_, m); ruu <- rep(NA_real_, m)
  if (n > 1L) {
    for (j in seq_len(m)) {
      rll[j] <- suppressWarnings(stats::cor(xa[, j, 1L], ya[, j, 1L]))
      ruu[j] <- suppressWarnings(stats::cor(xa[, j, 2L], ya[, j, 2L]))
    }
  } else {
    warning("correlation is not defined for a sample of size 1")
  }
  scalar <- 0.5 * trapz_alpha(rll + ruu, alpha)
  structure(list(alpha = alpha, r_lower = rll, r_upper = ruu, scalar = scalar),
            class = "FuzzyCorrelation")
}

#' Scalar covariance (Feng--Hu--Shu integral)
#'
#' \code{spectral_covariance} returns the scalar covariance
#' \deqn{\mathrm{Cov}_{FHS}(X,Y)=\tfrac12\int_0^1\bigl[\mathrm{Cov}(X^-(a),Y^-(a))+
#'        \mathrm{Cov}(X^+(a),Y^+(a))\bigr]\,da.}
#'
#' @param x,y Two \code{FuzzySample} objects (or \code{x} alone).
#' @param alpha Optional alpha grid.
#' @param ... Unused.
#' @return A length-one numeric.
#' @examples
#' set.seed(6)
#' x <- FuzzySample(replicate(30, triangular_fuzzy(0, 1, 2), simplify = FALSE))
#' spectral_covariance(x)
#' @export
spectral_covariance <- function(x, y = NULL, alpha = NULL, ...) {
  UseMethod("spectral_covariance")
}

#' @method spectral_covariance FuzzySample
#' @rdname spectral_covariance
#' @export
spectral_covariance.FuzzySample <- function(x, y = NULL, alpha = NULL, ...) {
  covariance(x, y, alpha)$scalar
}

#' Scalar variance (Feng--Hu--Shu integral)
#'
#' \code{spectral_variance} returns \code{spectral_covariance(x, x)}.
#'
#' @param x A \code{FuzzySample}.
#' @param alpha Optional alpha grid.
#' @param ... Unused.
#' @return A length-one numeric.
#' @examples
#' set.seed(7)
#' x <- FuzzySample(replicate(30, triangular_fuzzy(0, 1, 2), simplify = FALSE))
#' spectral_variance(x)
#' @export
spectral_variance <- function(x, alpha = NULL, ...) {
  spectral_covariance(x, x, alpha, ...)
}

#### FuzzyCovariance methods ####

#' @method print FuzzyCovariance
#' @export
print.FuzzyCovariance <- function(x, ...) {
  cat("<FuzzyCovariance>\n")
  cat(sprintf("  type    : %s\n", x$type))
  cat(sprintf("  levels  : %d\n", length(x$alpha)))
  cat(sprintf("  scalar  : %g\n", x$scalar))
  cat("  at levels alpha = 0, 0.5, 1:\n")
  idx <- round(stats::quantile(seq_along(x$alpha), c(0, 0.5, 1)))
  m <- matrix(NA_real_, length(idx), 4L,
              dimnames = list(format(x$alpha[idx]),
                              c("Var(lo)", "Cov(lo,up)", "Cov(up,lo)", "Var(up)")))
  for (k in seq_along(idx)) {
    cv <- x$cov[idx[k], , ]
    m[k, ] <- c(cv[1, 1], cv[1, 2], cv[2, 1], cv[2, 2])
  }
  print.default(m)
  invisible(x)
}

#' @method summary FuzzyCovariance
#' @export
summary.FuzzyCovariance <- function(object, ...) {
  a <- object$alpha
  df <- data.frame(alpha = a,
                   var_lower = object$cov[, 1L, 1L],
                   var_upper = object$cov[, 2L, 2L],
                   cov_cross = object$cov[, 1L, 2L],
                   scalar = object$scalar)
  structure(df, class = c("summary_FuzzyCovariance", "data.frame"), type = object$type)
}

#' @method print summary_FuzzyCovariance
#' @export
print.summary_FuzzyCovariance <- function(x, ...) {
  cat("<summary FuzzyCovariance>\n")
  cat(sprintf("  type   : %s\n  scalar: %g\n", attr(x, "type") %||% "?",
              x$scalar[1L]))
  print.data.frame(as.data.frame(x)[, c("alpha", "var_lower", "var_upper",
                                        "cov_cross")], row.names = FALSE)
  invisible(x)
}

#' @method as.data.frame FuzzyCovariance
#' @export
as.data.frame.FuzzyCovariance <- function(x, ...) {
  data.frame(alpha = x$alpha,
             var_lower = x$cov[, 1L, 1L],
             var_upper = x$cov[, 2L, 2L],
             cov_cross = x$cov[, 1L, 2L])
}

#' @method plot FuzzyCovariance
#' @export
plot.FuzzyCovariance <- function(x, main = "Level-wise covariance",
                                 which = c("variances", "ellipse", "all"),
                                 level = 0.5, ...) {
  which <- match.arg(which)
  a <- x$alpha
  if (which == "variances" || which == "all") {
    if (which == "all") {
      graphics::par(mfrow = c(2, 2))
      on.exit(graphics::par(mfrow = c(1, 1)), add = TRUE)
    }
    graphics::matplot(a, cbind(x$var_lower, x$var_upper, x$cov_cross),
                      type = "l", lty = 1:3, lwd = 2,
                      col = c("black", "grey40", "red"),
                      xlab = "alpha", ylab = "covariance",
                      main = paste(main, "(variances / cross)"))
    graphics::legend("topright",
                     legend = c("Var(lower)", "Var(upper)", "Cov(lo,up)"),
                     lty = 1:3, col = c("black", "grey40", "red"), bty = "n")
    if (which != "all") return(invisible(x))
  }
  if (which == "ellipse" || which == "all") {
    j <- which.min(abs(a - level))
    Sigma <- x$cov[j, , ]
    mu <- c(x$var_lower[j], x$var_upper[j])
    ell <- ellipse_points(Sigma, mu, level = 0.95)
    graphics::plot(ell[, 1L], ell[, 2L], type = "l", asp = 1,
                   xlab = "lower endpoint", ylab = "upper endpoint",
                   main = paste0("95% covariance ellipse (alpha = ", a[j], ")"))
    graphics::points(mu[1L], mu[2L], pch = 4, col = "red")
  }
  invisible(x)
}

#' Internal 2D ellipse points
#'
#' Internal helper used by plot.FuzzyCovariance.
#' @keywords internal
ellipse_points <- function(Sigma, mu, level = 0.95) {
  Sigma <- repair_cov_2x2(Sigma)
  ev <- eigen(Sigma, symmetric = TRUE)
  ang <- seq(0, 2 * pi, length.out = 200L)
  r <- sqrt(stats::qchisq(level, 2L))


  pts <- t(ev$vectors %*% diag(sqrt(ev$values)) %*%
           rbind(cos(ang) * r, sin(ang) * r))
  sweep(pts, 2L, mu, "+")
}

#### FuzzyCorrelation methods ####

#' @method print FuzzyCorrelation
#' @export
print.FuzzyCorrelation <- function(x, ...) {
  cat("<FuzzyCorrelation>\n")
  cat(sprintf("  levels : %d\n", length(x$alpha)))
  cat(sprintf("  scalar : %g\n", x$scalar))
  cat("  at levels alpha = 0, 0.5, 1:\n")
  idx <- round(stats::quantile(seq_along(x$alpha), c(0, 0.5, 1)))
  m <- rbind(lower = x$r_lower[idx], upper = x$r_upper[idx])
  print.default(m)
  invisible(x)
}

#' @method summary FuzzyCorrelation
#' @export
summary.FuzzyCorrelation <- function(object, ...) {
  df <- data.frame(alpha = object$alpha,
                   cor_lower = object$r_lower,
                   cor_upper = object$r_upper)
  structure(df, class = c("summary_FuzzyCorrelation", "data.frame"))
}

#' @method print summary_FuzzyCorrelation
#' @export
print.summary_FuzzyCorrelation <- function(x, ...) {
  cat("<summary FuzzyCorrelation>\n")
  print.data.frame(as.data.frame(x), row.names = FALSE)
  invisible(x)
}

#' @method as.data.frame FuzzyCorrelation
#' @export
as.data.frame.FuzzyCorrelation <- function(x, ...) {
  data.frame(alpha = x$alpha, cor_lower = x$r_lower, cor_upper = x$r_upper)
}

#' @method plot FuzzyCorrelation
#' @export
plot.FuzzyCorrelation <- function(x, main = "Level-wise correlation", ...) {
  graphics::matplot(x$alpha, cbind(x$r_lower, x$r_upper), type = "l",
                    lty = 1, lwd = 2, col = c("steelblue", "firebrick"),
                    xlab = "alpha", ylab = "correlation", ylim = c(-1, 1),
                    main = main, ...)
  graphics::abline(h = 0, lty = 3)
  graphics::legend("topright",
                   legend = c("cor(lower)", "cor(upper)"),
                   lty = 1, col = c("steelblue", "firebrick"), bty = "n")
  invisible(x)
}


