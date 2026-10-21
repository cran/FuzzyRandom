#### Simulation from the fuzzy normal distribution ####

#' Sample fuzzy random variables from a fuzzy normal distribution
#'
#' Draws \code{n} fuzzy numbers whose alpha-cut endpoints have the
#' bivariate-normal level-wise distribution \code{N_2(mu_alpha, Sigma_alpha)}.
#'
#' When \code{valid = TRUE} (the default) every draw is projected onto the
#' cone of valid, nested alpha-cuts (lower <= upper; lower non-decreasing and
#' upper non-increasing in \code{alpha}) so that it is a proper fuzzy number.
#' When \code{valid = FALSE} the raw, unprojected level-wise bivariate-normal
#' cuts are returned -- useful for checking the level-wise central limit
#' theorem, where projection would bias the marginals.
#'
#' @param n Sample size.
#' @param mean A \code{FuzzyNumber} (the fuzzy mean) or a \code{FuzzyNormal}.
#' @param cov A \code{2 x 2} covariance matrix, a length-4 vector, or a
#'   \code{FuzzyCovariance}. Ignored if \code{mean} is a \code{FuzzyNormal}.
#' @param valid Logical; whether to project draws to valid fuzzy numbers.
#' @param seed Optional integer seed (\code{NULL} leaves the RNG alone).
#' @return A \code{FuzzySample} of \code{n} fuzzy random variables.
#' @examples
#' set.seed(123)
#' s <- rnorm_fuzzy(80, triangular_fuzzy(0, 1, 2), diag(2) * 0.25)
#' plot(s, main = "fuzzy normal sample")
#' @export
rnorm_fuzzy <- function(n, mean, cov = NULL, valid = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (inherits(mean, "FuzzyNormal")) {
    fn <- mean
  } else {
    if (!inherits(mean, "FuzzyNumber"))
      stop("'mean' must be a FuzzyNumber or a FuzzyNormal")
    if (is.null(cov)) stop("please supply 'cov'")
    fn <- fuzzy_normal(mean, cov)
  }
  n <- as.integer(n)
  if (n < 1L) stop("'n' must be positive")
  alpha <- fn$alpha
  acov <- fn$cov$cov
  m <- length(alpha)
  arr <- array(NA_real_, c(n, m, 2L),
               dimnames = list(NULL, NULL, c("lower", "upper")))
  mu <- rbind(fn$mean$lower, fn$mean$upper)  # 2 x m
  for (j in seq_len(m)) {
    arr[, j, ] <- mvrnorm_2d(n, mu[, j, drop = TRUE], acov[j, , ])
  }
  if (isTRUE(valid)) {
    for (i in seq_len(n)) {
      pr <- project_to_valid_cuts(arr[i, , 1L], arr[i, , 2L], alpha)
      arr[i, , 1L] <- pr$lower
      arr[i, , 2L] <- pr$upper
    }
  }
  new_FuzzySample(arr, alpha)
}

#' @method simulate FuzzyNormal
#' @param object A \code{FuzzyNormal} (or \code{FuzzyCLT}) object.
#' @param nsim Number of simulations.
#' @param ... Unused.
#' @rdname rnorm_fuzzy
#' @export
simulate.FuzzyNormal <- function(object, nsim = 1, ...) {
  rnorm_fuzzy(nsim, object$mean, object$cov, valid = TRUE)
}

#' @method simulate FuzzyCLT
#' @rdname rnorm_fuzzy
#' @export
simulate.FuzzyCLT <- function(object, nsim = 1, ...) {
  ## draws of the sample mean, whose level cuts are N_2(mu, Sigma/n)
  rnorm_fuzzy(nsim, object$mean, object$cov, valid = TRUE)
}
