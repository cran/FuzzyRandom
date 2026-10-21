test_that("triangular_fuzzy constructs a valid number", {
  fn <- triangular_fuzzy(0, 1, 2)
  expect_s3_class(fn, "FuzzyNumber")
  expect_true(all(fn$lower <= fn$upper))
  expect_true(all(diff(fn$lower) >= 0))
  expect_true(all(diff(fn$upper) <= 0))
  expect_equal(fn$lower[1], 0)
  expect_equal(fn$upper[1], 2)
  expect_equal(fn$lower[length(fn$alpha)], 1)
  expect_equal(fn$upper[length(fn$alpha)], 1)
})

test_that("FuzzySample stores valid cuts", {
  set.seed(1)
  fs <- FuzzySample(replicate(25, triangular_fuzzy(0, 1, 2), simplify = FALSE))
  expect_s3_class(fs, "FuzzySample")
  expect_equal(n_obs(fs), 25)
  expect_equal(dim(fs$x), c(25, length(fs$alpha), 2))
  expect_true(all(fs$x[, , 1] <= fs$x[, , 2]))
})

test_that("expectation and variance are level-wise", {
  set.seed(3)
  fs <- rnorm_fuzzy(80, triangular_fuzzy(0, 1, 2), diag(2) * 0.25)
  ex <- expectation(fs)
  expect_s3_class(ex, "FuzzyNumber")
  expect_true(all(ex$lower <= ex$upper))
  cv <- variance(fs)
  expect_s3_class(cv, "FuzzyCovariance")
  expect_true(all(cv$cov[, 1, 1] >= 0))
  expect_true(all(cv$cov[, 2, 2] >= 0))
  expect_gt(cv$scalar[1], 0)
})

test_that("covariance and correlation are well defined", {
  set.seed(5)
  a <- rnorm_fuzzy(40, triangular_fuzzy(0, 1, 2), diag(2) * 0.25)
  b <- rnorm_fuzzy(40, triangular_fuzzy(1, 2, 3), diag(2) * 0.25)
  cv <- covariance(a, b)
  expect_s3_class(cv, "FuzzyCovariance")
  cr <- correlation(a, b)
  expect_true(all(cr$r_lower >= -1 & cr$r_lower <= 1))
  expect_true(all(cr$r_upper >= -1 & cr$r_upper <= 1))
})

test_that("fuzzy_ecdf is monotonic and bounded", {
  set.seed(4)
  fs <- FuzzySample(replicate(40, triangular_fuzzy(0, 1, 2), simplify = FALSE))
  fe <- fuzzy_ecdf(fs, alpha = c(0, 0.5, 1))
  expect_s3_class(fe, "FuzzyECDF")
  expect_equal(nrow(fe$F_lower), 3)
  for (j in seq_len(nrow(fe$F_lower))) {
    expect_gte(min(fe$F_lower[j, ]), 0)
    expect_lte(max(fe$F_lower[j, ]), 1)
    expect_true(all(diff(fe$F_lower[j, ]) >= 0))
  }
})

test_that("rnorm_fuzzy returns a valid sample", {
  set.seed(12)
  s <- rnorm_fuzzy(60, triangular_fuzzy(0, 1, 2), diag(2) * 0.25)
  expect_s3_class(s, "FuzzySample")
  expect_equal(s$n, 60)
  expect_true(all(s$x[, , 1] <= s$x[, , 2]))
})
