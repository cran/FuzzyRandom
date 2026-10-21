## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 4, dpi = 80)
library(FuzzyRandom)
set.seed(2024)

## ----fuzzy-numbers------------------------------------------------------------
tri <- triangular_fuzzy(0, 1, 2)
trap <- trapezoidal_fuzzy(-1, 0, 1, 3)
plot(tri)

## ----sample-------------------------------------------------------------------
s <- FuzzySample(replicate(60, triangular_fuzzy(0, 1, 2), simplify = FALSE))
plot(s)

## ----moments------------------------------------------------------------------
expectation(s)
variance(s)
correlation(s, FuzzySample(replicate(60, triangular_fuzzy(1, 2, 3),
                                  simplify = FALSE)))

## ----ecdf-ks------------------------------------------------------------------
fe <- fuzzy_ecdf(s, alpha = c(0, 0.5, 1))
plot(fe)
ks_test(s, p0 = "pnorm")

## ----clt----------------------------------------------------------------------
cl <- clt(s)
summary(cl)
plot(simulate(cl, nsim = 80))

## ----lln----------------------------------------------------------------------
w <- lln(s)
plot(w)

