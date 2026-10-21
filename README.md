# FuzzyRandom

Level-wise (*Puri--Ralescu*) statistics for fuzzy random variables.

Fuzzy numbers are represented through their alpha-level cuts
`[X^-(alpha), X^+(alpha)]`, and every statistical operator acts **level-wise**
so that the expectation, variance, covariance, empirical distribution,
goodness-of-fit, central limit theorem and law of large numbers all live in a
single coherent framework.

## Installation

Install from source (after `R CMD build` / `R CMD INSTALL`):

```r
# install.packages("remotes")
remotes::install_dir("path/to/FuzzyRandom")
```

## Quick start

```r
library(FuzzyRandom)
set.seed(123)
s <- rnorm_fuzzy(80, triangular_fuzzy(0, 1, 2), diag(2) * 0.25)

plot(s)
plot(fuzzy_ecdf(s))
ks_test(s, p0 = "pnorm")

cl <- clt(s)
summary(cl)
plot(s <- simulate(cl, nsim = 60))
```
