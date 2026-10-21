#### Fuzzy empirical distribution function ####
####
#### For a sample X_1,...,X_n of fuzzy random variables and a level alpha the
#### empirical distribution of the alpha-cut is the pair of empirical CDFs of
#### its endpoints:
####   F_hat^-_alpha(x) = (1/n) sum_i 1_{X_i^-(alpha) <= x}
####   F_hat^+_alpha(x) = (1/n) sum_i 1_{X_i^+(alpha) <= x}
#### (Hesamian & Taheri, Kybernetika 49(6):962-982, 2013). By the classical
#### Glivenko-Cantelli theorem applied to each endpoint sequence, F_hat^-_alpha
#### and F_hat^+_alpha converge to the population level-wise CDFs at every
#### alpha (Puri & Ralescu, 1986).

new_FuzzyECDF <- function(alpha, at, F_lower, F_upper, n) {
  structure(list(alpha = alpha, x = at,
                 F_lower = F_lower, F_upper = F_upper, n = n),
            class = "FuzzyECDF")
}

#' Empirical distribution function of a fuzzy random sample
#'
#' Level-wise empirical cumulative distribution functions of the alpha-cut
#' endpoints.
#'
#' @param x A \code{FuzzySample}.
#' @param alpha Numeric vector of levels (default: the levels stored in
#'   \code{x}).
#' @param at Numeric grid on which to evaluate the CDFs (default: 201 points
#'   over the data range).
#' @param ... Unused.
#' @return An object of class \code{"FuzzyECDF"}.
#' @examples
#' set.seed(8)
#' s <- FuzzySample(replicate(40, triangular_fuzzy(0, 1, 2), simplify = FALSE))
#' fe <- fuzzy_ecdf(s, alpha = c(0, 0.25, 0.5, 0.75, 1))
#' plot(fe)
#' @export
fuzzy_ecdf <- function(x, alpha = NULL, at = NULL, ...) {
  if (is.null(alpha)) alpha <- x$alpha
  a <- standardise_alpha(alpha)
  UseMethod("fuzzy_ecdf")
}

#' @method fuzzy_ecdf FuzzySample
#' @rdname fuzzy_ecdf
#' @export
fuzzy_ecdf.FuzzySample <- function(x, alpha = NULL, at = NULL, ...) {
  if (is.null(alpha)) a <- x$alpha else a <- standardise_alpha(alpha)
  arr <- interpolate_sample(x, a)
  m <- length(a); n <- x$n
  if (is.null(at)) {
    rng <- range(c(arr[, , 1L], arr[, , 2L]))
    at <- seq(rng[1L], rng[2L], length.out = 201L)
  }
  at <- as.numeric(at)
  Fl <- matrix(NA_real_, m, length(at))
  Fu <- matrix(NA_real_, m, length(at))
  for (j in seq_len(m)) {
    fLj <- stats::ecdf(arr[, j, 1L])
    fUj <- stats::ecdf(arr[, j, 2L])
    Fl[j, ] <- fLj(at)
    Fu[j, ] <- fUj(at)
  }
  new_FuzzyECDF(a, at, Fl, Fu, n)
}

#' @method print FuzzyECDF
#' @export
print.FuzzyECDF <- function(x, ...) {
  cat("<FuzzyECDF>\n")
  cat(sprintf("  n      : %d\n", x$n))
  cat(sprintf("  levels : %d\n", length(x$alpha)))
  cat(sprintf("  grid   : %d points\n", length(x$x)))
  invisible(x)
}

#' @method summary FuzzyECDF
#' @export
summary.FuzzyECDF <- function(object, ...) {
  a <- object$alpha
  Fl <- object$F_lower
  Fu <- object$F_upper
  at <- object$x
  lo <- sapply(seq_along(a), function(j)
    at[which(Fl[j, ] > 0)[1L]])
  hi <- sapply(seq_along(a), function(j)
    at[which(Fu[j, ] < 1)[1L]])
  df <- data.frame(alpha = a, support_lo = lo, support_hi = hi)
  structure(df, class = c("summary_FuzzyECDF", "data.frame"))
}

#' @method print summary_FuzzyECDF
#' @export
print.summary_FuzzyECDF <- function(x, ...) {
  cat("<summary FuzzyECDF>\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}

#' @method as.data.frame FuzzyECDF
#' @export
as.data.frame.FuzzyECDF <- function(x, ...) {
  m <- length(x$alpha); k <- length(x$x)
  do.call(rbind, lapply(seq_len(m), function(j)
    data.frame(alpha = x$alpha[j], x = x$x,
               F_lower = x$F_lower[j, ], F_upper = x$F_upper[j, ])))
}

#' @method [ FuzzyECDF
#' @export
'[.FuzzyECDF' <- function(x, i, ...) {
  new_FuzzyECDF(x$alpha[i], x$x,
                x$F_lower[i, , drop = FALSE], x$F_upper[i, , drop = FALSE], x$n)
}

#' @method plot FuzzyECDF
#' @export
plot.FuzzyECDF <- function(x, main = "Fuzzy empirical CDF",
                           levels = NULL, ...) {
  a <- x$alpha; k <- length(a)
  if (is.null(levels)) levels <- sort(unique(c(0, 0.25, 0.5, 0.75, 1)))
  idx <- sapply(levels, function(v) which.min(abs(a - v)))
  idx <- sort(unique(pmin(pmax(idx, 1L), k)))
  cols <- grDevices::rainbow(length(idx))
  graphics::plot(range(x$x), c(0, 1), type = "n", xlab = "x", ylab = "F(x)",
                 main = main, ...)
  for (q in seq_along(idx)) {
    j <- idx[q]
    graphics::lines(x$x, x$F_lower[j, ], col = cols[q], lwd = 2)
    graphics::lines(x$x, x$F_upper[j, ], col = cols[q], lty = 2, lwd = 2)
  }
  graphics::legend("topleft", legend = paste0("alpha = ", a[idx]),
                   col = cols, lty = 1, lwd = 2, bty = "n")
  graphics::abline(h = c(0, 1), lty = 3)
  invisible(x)
}

#### Level-wise Kolmogorov-Smirnov tests ####

#' One- and two-sample fuzzy Kolmogorov-Smirnov tests
#'
#' Performs, at each level \code{alpha}, a Kolmogorov-Smirnov test on the
#' lower and upper endpoints of the alpha-cuts. With \code{p0} a one-sample
#' test against the CDF \code{p0} is run; otherwise a two-sample test against
#' the second sample \code{y} is run. The result is a family of level-wise
#' tests consistent with the Puri-Ralescu level-wise approach.
#'
#' @param x A \code{FuzzySample}.
#' @param y A second \code{FuzzySample} for a two-sample test (default NULL).
#' @param p0 Null CDF for the one-sample test: a character string naming a CDF
#'   in \code{stats} (e.g. \code{"pnorm"}) or a function \code{p0(q)}.
#' @param alpha Numeric vector of levels (default: those of \code{x}).
#' @param ... Additional arguments passed on to \code{\link[stats]{ks.test}}.
#' @return An object of class \code{"FuzzyKS"}.
#' @examples
#' set.seed(9)
#' s <- FuzzySample(replicate(50, triangular_fuzzy(-1, 0, 1), simplify = FALSE))
#' ks_test(s, p0 = "pnorm")
#' @export
ks_test <- function(x, ...) UseMethod("ks_test")

new_FuzzyKS <- function(alpha, n, D_lower, D_upper, p_lower, p_upper,
                        method, data.name) {
  structure(
    list(alpha = alpha, n = n,
           D_lower = D_lower, D_upper = D_upper,
           p_lower = p_lower, p_upper = p_upper,
           method = method, data.name = data.name),
    class = "FuzzyKS")
}

#' @method ks_test FuzzySample
#' @rdname ks_test
#' @export
ks_test.FuzzySample <- function(x, y = NULL, p0 = NULL, alpha = NULL, ...) {
  if (is.null(alpha)) a <- x$alpha else a <- standardise_alpha(alpha)
  xa <- interpolate_sample(x, a)
  m <- length(a); n <- x$n
  Dl <- numeric(m); Du <- numeric(m); pl <- numeric(m); pu <- numeric(m)
  if (is.null(y) && !is.null(p0)) {
    method <- paste0("Level-wise one-sample KS test (p0 = ",
                     if (is.character(p0)) p0 else "function", ")")
    for (j in seq_len(m)) {
      t1 <- suppressWarnings(stats::ks.test(xa[, j, 1L], p0, ...))
      t2 <- suppressWarnings(stats::ks.test(xa[, j, 2L], p0, ...))
      Dl[j] <- t1$statistic; Du[j] <- t2$statistic
      pl[j] <- t1$p.value;  pu[j] <- t2$p.value
    }
  } else if (!is.null(y)) {
    if (!inherits(y, "FuzzySample")) stop("'y' must be a FuzzySample")
    ya <- interpolate_sample(y, a)
    method <- "Level-wise two-sample KS test"
    for (j in seq_len(m)) {
      t1 <- suppressWarnings(stats::ks.test(xa[, j, 1L], ya[, j, 1L], ...))
      t2 <- suppressWarnings(stats::ks.test(xa[, j, 2L], ya[, j, 2L], ...))
      Dl[j] <- t1$statistic; Du[j] <- t2$statistic
      pl[j] <- t1$p.value;  pu[j] <- t2$p.value
    }
  } else {
    stop("either 'p0' (one-sample) or 'y' (two-sample) must be supplied")
  }
  new_FuzzyKS(a, n, Dl, Du, pl, pu, method,
              paste(deparse(substitute(x)),
                    if (!is.null(y)) deparse(substitute(y))))
}

#' @method print FuzzyKS
#' @export
print.FuzzyKS <- function(x, ...) {
  cat("<FuzzyKS>\n")
  cat(x$method, "\n")
  cat(sprintf("  levels: %d\n", length(x$alpha)))
  m <- length(x$alpha)
  idx <- round(stats::quantile(seq_len(m), c(0, 0.5, 1)))
  idx <- sort(unique(pmax(1L, pmin(m, idx))))
  Dsel <- pmax(x$D_lower[idx], x$D_upper[idx])
  cat("  D (max(lower, upper)) at selected levels:",
      paste(format(Dsel), collapse = " "), "\n")
  cat(sprintf("  scalar p-value (min over levels): %g\n",
              min(x$p_lower, x$p_upper)))
  invisible(x)
}

#' @method summary FuzzyKS
#' @export
summary.FuzzyKS <- function(object, ...) {
  df <- data.frame(alpha = object$alpha,
                   D_lower = object$D_lower, D_upper = object$D_upper,
                   p_lower = object$p_lower, p_upper = object$p_upper)
  structure(df, class = c("summary_FuzzyKS", "data.frame"))
}

#' @method print summary_FuzzyKS
#' @export
print.summary_FuzzyKS <- function(x, ...) {
  cat("<summary FuzzyKS>\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}

#' @method as.data.frame FuzzyKS
#' @export
as.data.frame.FuzzyKS <- function(x, ...) {
  data.frame(alpha = x$alpha, D_lower = x$D_lower, D_upper = x$D_upper,
             p_lower = x$p_lower, p_upper = x$p_upper)
}

#' @method plot FuzzyKS
#' @export
plot.FuzzyKS <- function(x, main = "KS statistic vs level", ...) {
  D <- pmax(x$D_lower, x$D_upper)
  graphics::plot(x$alpha, D, type = "b", pch = 19, lwd = 2,
                 xlab = "alpha", ylab = "D", ylim = c(0, 1), main = main, ...)
  graphics::abline(h = 0.05, lty = 2, col = "red")
  invisible(x)
}
