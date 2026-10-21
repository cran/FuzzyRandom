#### Package-level documentation ####

#' Statistics of Fuzzy Random Variables: a Puri-Ralescu Framework
#'
#' `FuzzyRandom` is a computationally self-contained implementation of the
#' level-wise ("Puri-Ralescu") theory of statistics for fuzzy random
#' variables.
#'
#' Fuzzy numbers are represented through their `alpha`-level sets
#' `[X^-(alpha), X^+(alpha)]`, and every statistical operator is defined
#' level-wise, so that the **expectation**, **variance**, **covariance**,
#' **empirical distribution**, **law of large numbers** and **central limit
#' theorem** are all handled in a single coherent framework.
#'
#' The key design choices are:
#'
#' * A *fuzzy number* is a list of `alpha` levels and the matching
#'   non-decreasing `lower` / non-increasing `upper` envelopes
#'   (class `FuzzyNumber`).
#' * A *sample of fuzzy random variables* is an `(n, m, 2)` array of cuts
#'   (class `FuzzySample`), `m` being the number of levels.
#' * The covariance of two fuzzy random variables is returned as a
#'   level-wise `2 x 2` matrix of covariances between the lower and upper
#'   endpoints (class `FuzzyCovariance`); scalar reductions are provided
#'   through the Feng-Hu-Shu integral and scalar product.
#' * The central limit theorem yields the *fuzzy normal* distribution
#'   (class `FuzzyNormal`), whose level cuts are bivariate-normal intervals.
#'
#' The package only depends on base R (`stats`, `graphics`).
#'
#' @name FuzzyRandom
#' @import stats graphics grDevices
#' @references
#' Puri, M. L. and Ralescu, D. A. (1986). Fuzzy random variables.
#'   \emph{Journal of Mathematical Analysis and Applications},
#'   \strong{114}(2), 409-422.
#'
#' Klement, E. P., Puri, M. L. and Ralescu, D. A. (1986). Limit theorems for
#'   fuzzy random variables. \emph{Proceedings of the Royal Society A},
#'   \strong{407}, 171-182.
#'
#' Feng, Y., Hu, L. and Shu, H. (2002). The variance and covariance of fuzzy
#'   random variables and their applications. \emph{Fuzzy Sets and Systems},
#'   \strong{128}, 335-345.
"_PACKAGE"
