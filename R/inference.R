#### Fuzzy normal distribution, CLT and LLN (Puri-Ralescu framework) ####
####
#### A fuzzy normal is specified by a fuzzy mean (FuzzyNumber) and a level-wise
#### 2x2 covariance matrix Sigma(alpha) = [[var^-(alpha), cov(alpha)],
#### [cov(alpha), var^+(alpha)]]. At each level the alpha-cut endpoints are
#### jointly normal with this covariance. By the level-wise multivariate CLT
#### (Puri & Ralescu, 1986), the sample-mean fuzzy number satisfies
####   sqrt(n) * (Xbar - mu)(alpha)  -->  N_2(0, Sigma(alpha))
#### so the covariance of Xbar is Sigma(alpha)/n at level alpha.

#' Construct a fuzzy normal distribution
#'
#' A fuzzy normal distribution is the family, over levels \code{alpha}, of
#' bivariate-normal laws on the alpha-cut endpoints with fuzzy mean \code{mean}
#' and level-wise covariance \code{cov}.
#'
#' @param mean A \code{FuzzyNumber} (the fuzzy mean).
#' @param cov A \code{2 x 2} covariance matrix, a length-4 vector, or a
#'   \code{FuzzyCovariance} whose levels match \code{mean}.
#' @param n Optional sample size associated with the mean (used by
#'   \code{\link{clt}}).
#' @return An object of class \code{"FuzzyNormal"}.
#' @examples
#' fn <- fuzzy_normal(triangular_fuzzy(0, 1, 2), diag(2))
#' simulate(fn, nsim = 5)
#' @export
fuzzy_normal <- function(mean, cov = NULL, n = NULL) {
  if (inherits(mean, "FuzzyNormal")) return(mean)
  if (!inherits(mean, "FuzzyNumber"))
    stop("'mean' must be a FuzzyNumber")
  alpha <- mean$alpha
  m <- length(alpha)
  if (inherits(cov, "FuzzyCovariance")) {
    if (length(cov$alpha) != m)
      stop("alpha grid of 'cov' does not match 'mean'")
    cv_arr <- cov$cov
    cv_scalar <- cov$scalar
  } else {
    if (is.null(cov)) stop("please supply 'cov'")
    cm <- if (is.matrix(cov)) cov else matrix(as.numeric(cov), 2L, 2L, byrow = TRUE)
    if (nrow(cm) != 2L || ncol(cm) != 2L)
      stop("'cov' must be a 2 x 2 matrix (or length-4 vector)")
    cv_arr <- array(rep(as.vector(cm), m), c(2L, 2L, m))
    cv_arr <- aperm(cv_arr, c(3L, 1L, 2L))
    cv_scalar <- trapz_alpha(0.5 * (cm[1L, 1L] + cm[2L, 2L]), alpha)
    cv_scalar <- rep(cv_scalar, m)
  }
  fc <- new_FuzzyCovariance(alpha, cv_arr, cv_scalar, "covariance")
  structure(list(mean = mean, cov = fc, sigma = fc, n = n, alpha = alpha),
            class = "FuzzyNormal")
}

#' @method print FuzzyNormal
#' @export
print.FuzzyNormal <- function(x, ...) {
  cls <- if ("FuzzyCLT" %in% class(x)) "Fuzzy CLT" else "Fuzzy normal"
  cat("<", cls, ">\n", sep = "")
  cat(sprintf("  mean support : [%g, %g]\n",
              min(x$mean$lower), max(x$mean$upper)))
  cat(sprintf("  levels       : %d\n", length(x$alpha)))
  cat(sprintf("  scalar var   : %g\n", mean(x$cov$scalar)))
  if (!is.null(x$n)) cat(sprintf("  n            : %d\n", x$n))
  cat(sprintf("  method       : %s\n",
              if ("FuzzyCLT" %in% class(x)) x$method else "fuzzy_normal"))
  invisible(x)
}

#' @method summary FuzzyNormal
#' @export
summary.FuzzyNormal <- function(object, ...) {
  a <- object$alpha
  S <- object$cov$cov
  df <- data.frame(alpha = a,
                   mean_lower = object$mean$lower,
                   mean_upper = object$mean$upper,
                   var_lower = S[, 1L, 1L],
                   var_upper = S[, 2L, 2L],
                   cov_cross = S[, 1L, 2L],
                   scalar = object$cov$scalar)
  structure(df, class = c("summary_FuzzyNormal", "data.frame"))
}

#' @method print summary_FuzzyNormal
#' @export
print.summary_FuzzyNormal <- function(x, ...) {
  cat("<summary FuzzyNormal>\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}

#' @method as.data.frame FuzzyNormal
#' @export
as.data.frame.FuzzyNormal <- function(x, ...) as.data.frame(summary(x))

#' @method expectation FuzzyNormal
#' @export
expectation.FuzzyNormal <- function(x, alpha = NULL, ...) x$mean

#' @method variance FuzzyNormal
#' @export
variance.FuzzyNormal <- function(x, alpha = NULL, ...) x$cov

#' @method plot FuzzyNormal
#' @export
plot.FuzzyNormal <- function(x, k = 2, main = "Fuzzy normal", ...) {
  a <- x$alpha
  mu <- x$mean
  S <- x$cov$cov
  se_lo <- sqrt(pmax(S[, 1L, 1L], 0))
  se_up <- sqrt(pmax(S[, 2L, 2L], 0))
  lo <- mu$lower - k * se_lo
  up <- mu$upper + k * se_up
  pr <- project_to_valid_cuts(lo, up, a)
  plot(mu, main = main, ...)
  graphics::polygon(c(pr$lower, rev(pr$upper)), c(a, rev(a)),
                    col = grDevices::rgb(0.92, 0.40, 0.15, 0.22), border = NA)
  graphics::lines(pr$lower, a, col = grDevices::rgb(0.86, 0.15, 0.16), lty = 2)
  graphics::lines(pr$upper, a, col = grDevices::rgb(0.86, 0.15, 0.16), lty = 2)
  invisible(x)
}

#### Level-wise central limit theorem ####

#' Central limit theorem for a fuzzy random sample
#'
#' Uses the level-wise central limit theorem (Puri & Ralescu, 1986): at each
#' level \code{alpha} the sample-mean endpoints are approximately normal, so
#' the fuzzy sample mean \code{Xbar} is approximately a
#' \code{\link{fuzzy_normal}} with mean \code{\link{expectation}(x)} and level-
#' wise covariance \code{\link{variance}(x) / n}.
#'
#' @param x A \code{FuzzySample}.
#' @param alpha Numeric levels (default: those of \code{x}).
#' @param ... Unused.
#' @return An object of class \code{c("FuzzyCLT","FuzzyNormal")}. Its
#'   \code{cov} is the covariance of \code{Xbar} (i.e. \code{Sigma / n});
#'   \code{sigma} is the limiting covariance \code{Sigma}.
#' @examples
#' set.seed(7)
#' s <- rnorm_fuzzy(500, triangular_fuzzy(0, 1, 2), diag(2) * 0.25)
#' cl <- clt(s); summary(cl)
#' simulate(cl, nsim = 4)
#' @export
clt <- function(x, alpha = NULL, ...) {
  if (!inherits(x, "FuzzySample"))
    stop("'x' must be a FuzzySample")
  if (is.null(alpha)) a <- x$alpha else a <- standardise_alpha(alpha)
  mu <- expectation(x, a)            # FuzzyNumber
  sigma <- variance(x, a)           # FuzzyCovariance (Sigma)
  n <- x$n
  cv <- sigma$cov / n
  fc <- new_FuzzyCovariance(sigma$alpha, cv, sigma$scalar / n, "covariance")
  structure(list(mean = mu, cov = fc, sigma = sigma, n = n,
                 alpha = sigma$alpha, method = "Central limit theorem"),
            class = c("FuzzyCLT", "FuzzyNormal"))
}

#' @method print FuzzyCLT
#' @export
print.FuzzyCLT <- function(x, ...) {
  cat("<Fuzzy CLT>\n")
  cat(sprintf("  mean support : [%g, %g]\n",
              min(x$mean$lower), max(x$mean$upper)))
  cat(sprintf("  levels       : %d\n", length(x$alpha)))
  cat(sprintf("  n            : %d\n", x$n))
  cat(sprintf("  scalar var   : %g\n", mean(x$cov$scalar)))
  cat(sprintf("  scalar limit : %g\n", mean(x$sigma$scalar)))
  invisible(x)
}

#' @method summary FuzzyCLT
#' @export
summary.FuzzyCLT <- function(object, ...) {
  a <- object$alpha; S <- object$cov$cov
  df <- data.frame(alpha = a, mean_lower = object$mean$lower,
                   mean_upper = object$mean$upper, var_lower = S[, 1L, 1L],
                   var_upper = S[, 2L, 2L], cov_cross = S[, 1L, 2L],
                   scalar = object$cov$scalar,
                   scalar_limit = object$sigma$scalar)
  structure(df, class = c("summary_FuzzyCLT", "data.frame"))
}

#' @method print summary_FuzzyCLT
#' @export
print.summary_FuzzyCLT <- function(x, ...) {
  cat("<summary FuzzyCLT>\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}

#### Level-wise law of large numbers ####

#' Level-wise sample-mean convergence (LLN)
#'
#' Following Puri & Ralescu (1986), statistics on fuzzy random variables are
#' evaluated level-wise on the alpha-cuts. This function computes the fuzzy
#' sample mean for an increasing sequence of sample sizes and measures, at a
#' fixed level, the sup-norm distance between the size-\code{s} mean and the
#' full-sample mean. By the level-wise (fuzzy) law of large numbers this
#' distance shrinks to zero with the sample size (at rate \code{1/sqrt(s)}
#' under the level-wise CLT).
#'
#' @param x A \code{FuzzySample}.
#' @param sizes Integer vector of sample sizes (default: a grid from 5 up to
#'   \code{x$n}).
#' @param alpha Optional alpha grid on which to evaluate the means.
#' @param levels Numeric levels at which to measure convergence, used by the
#'   \code{summary}, \code{as.data.frame} and \code{plot} methods
#'   (default \code{c(0.5, 1)}).
#' @param ... Unused.
#' @return An object of class \code{"FuzzyLLN"} with components \code{sizes},
#'   \code{means}, \code{full} and \code{alpha}.
#' @examples
#' set.seed(23)
#' s <- rnorm_fuzzy(400, triangular_fuzzy(0, 1, 2), diag(2) * 0.25)
#' w <- lln(s)
#' plot(w)
#' summary(w)
#' @export
lln <- function(x, sizes = NULL, alpha = NULL, levels = c(0.5, 1), ...) {
  if (!inherits(x, "FuzzySample"))
    stop("lln() expects a FuzzySample")
  if (is.null(alpha)) alpha <- x$alpha else alpha <- standardise_alpha(alpha)
  n <- x$n
  if (is.null(sizes)) sizes <- round(seq(5, n, length.out = 8L))
  sizes <- as.integer(sort(unique(pmax(2L, pmin(as.integer(sizes), n)))))
  means <- lapply(sizes, function(s) expectation(x[seq_len(s)], alpha))
  full <- expectation(x, alpha)
  structure(list(sizes = sizes, means = means, full = full, alpha = alpha,
                 levels = levels, n = n),
            class = "FuzzyLLN")
}

#' @method print FuzzyLLN
#' @export
print.FuzzyLLN <- function(x, ...) {
  cat("<FuzzyLLN>\n")
  cat(sprintf("  n            : %d\n", x$n))
  cat(sprintf("  sizes used   : %d\n", length(x$sizes)))
  cat(sprintf("  size range   : %d .. %d\n", min(x$sizes), max(x$sizes)))
  cat("  levels       :", paste(x$levels, collapse = ", "), "\n", sep = " ")
  invisible(x)
}

#' Sup-norm distance of each mean to the full mean, at the requested levels
#'
#' @param means List of \code{FuzzyNumber} means (size-ascending).
#' @param full The full-sample \code{FuzzyNumber} mean.
#' @param alpha Alpha grid shared by \code{means} and \code{full}.
#' @param levels Levels at which the distance is measured.
#' @return A matrix with one column per level.
#' @keywords internal
lln_distances <- function(means, full, alpha, levels) {
  out <- matrix(NA_real_, length(means), length(levels))
  for (q in seq_along(levels)) {
    j <- which.min(abs(alpha - levels[q]))
    out[, q] <- vapply(means, function(m)
      max(abs(m$lower[j] - full$lower[j]),
          abs(m$upper[j] - full$upper[j])),
      numeric(1))
  }
  colnames(out) <- paste0("alpha=", levels)
  out
}

#' @method summary FuzzyLLN
#' @export
summary.FuzzyLLN <- function(object, ...) {
  d <- lln_distances(object$means, object$full, object$alpha, object$levels)
  df <- data.frame(size = object$sizes, d, check.names = FALSE)
  structure(df, class = c("summary_FuzzyLLN", "data.frame"))
}

#' @method print summary_FuzzyLLN
#' @export
print.summary_FuzzyLLN <- function(x, ...) {
  cat("<summary FuzzyLLN>\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}

#' @method as.data.frame FuzzyLLN
#' @export
as.data.frame.FuzzyLLN <- function(x, ...) {
  as.data.frame(summary(x), optional = TRUE)
}

#' @method plot FuzzyLLN
#' @export
plot.FuzzyLLN <- function(x, main = "LLN: mean distance vs sample size", ...) {
  levels <- x$levels
  cols <- grDevices::rainbow(length(levels))
  d <- lln_distances(x$means, x$full, x$alpha, levels)
  s <- x$sizes
  graphics::plot(s, d[, 1L], type = "n", log = "xy",
                 xlim = range(s), ylim = range(d, na.rm = TRUE),
                 xlab = "sample size", ylab = "distance of mean to full",
                 main = main, ...)
  for (q in seq_along(levels))
    graphics::lines(s, d[, q], col = cols[q], lwd = 2, type = "b", pch = 19)
  graphics::lines(s, d[1L, 1L] * sqrt(s[1L] / s), lty = 2, col = "grey40")
  graphics::legend("topright", legend = paste0("level ", levels),
                   col = cols, lty = 1, lwd = 2, pch = 19, bty = "n")
  invisible(x)
}
