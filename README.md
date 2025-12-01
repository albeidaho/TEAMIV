# TEAMIV
Required functions for TEAM-IV method (Team-based Evaluation and Aggregation with MCP for Instrumental Variables)

## Installation
source(
  "https://raw.githubusercontent.com/albeidaho/TEAMIV/9c7cba9/TEAMIV_functions.R"
)

## Dependencies

This script assumes the following R packages are installed and loaded:

- dplyr
- tibble
- purrr
- glmnet
- ncvreg

## Input of TEAM_IV() function
- `y`: numeric vector of outcomes (length `n`).
- `d`: numeric vector of exposures (length `n`).
- `Z`: numeric matrix of candidate instruments (`n x L`), with each column corresponding to one instrument.
- `K`: integer; number of folds for cross-validation used in both team construction and RMSE calculation (default: `20`).
- `lam_se`: nonnegative scalar; SE-multiplier in the “X–SE rule” used to choose λ for the final MNet/MCP fits (does **not** affect team construction; default: `0.1`).
- `check_gam_sign`: logical; if `TRUE`, TEAM-IV flips the sign of columns of `Z` whose ridge-estimated loading (`gamma_hat`) is negative before constructing `d_adj` (default: `FALSE`).
- `verbose`: integer; if `>= 1`, prints the teamness matrix and the results tibble to the console (default: `0`).
- `seed`: optional integer; if non-`NULL`, used to set the random seed (for reproducible folds and preprocessing).

## Output of TEAM_IV() function
`TEAM_IV()` returns a tibble with one row per team (and per union of teams), containing:

- `iv`: identifier for the instrument cluster (team) or union of clusters used in computation of beta-hat.
- `valid_ivs`: integer vector (list-column) giving the indices of instruments in `Z` (original ordering) treated as valid for this configuration.
- `n_valid`: number of valid instruments (`length(valid_ivs)`).
- `mean_rmse`: TEAM-IV/sisVIVE-style cross-validated root mean squared error (RMSE) for this configuration.
- `se_rmse`: standard error of `mean_rmse` across folds.
- `beta_afl`: beta-hat using the MNet estimator with mixing parameter chosen (between 0.5 and 0.75) to minimize cross-validated error (not generally recommended for reporting).
- `beta_a50`: beta-hat using MNet with mixing parameter 0.5.
- `beta_a75`: beta-hat using MNet with mixing parameter 0.75.
- `beta_mcp`: beta-hat using MCP (alpha = 1).
- `fold_mse`: numeric vector (list-column) of fold-wise mean squared errors from the sisVIVE-style cross-validation loss.


## Example usage:
```
library(MASS)
source(
  "https://raw.githubusercontent.com/albeidaho/TEAMIV/9c7cba9/TEAMIV_functions.R"
)
n = 2000 # sample size
L = 10 # number of candidate IVs
s = 6 # number of invalid IVs
alpha_star_const = c(0.6)
cov_offdiag = c(0.8)
beta_star = c(1L)
cov_mat =
  matrix(cov_offdiag, nrow = 2, ncol = 2) +
  diag(1 - cov_offdiag, nrow = 2, ncol = 2)
Z = MASS::mvrnorm(n, mu = rep(0, L), Sigma = diag(1, L))
# instruments
gamma_star = runif(L, 0.08, 0.14)
mu = Z %*% gamma_star # exposure means
errors = MASS::mvrnorm(n, mu = rep(0, 2), Sigma = cov_mat)
colnames(errors) = c('epsilon', 'nu')
d = mu + errors[, 'nu'] #exposure values
alpha_star =
  c(
    rep(1, s),
    rep(0, L - s)
  ) *
  alpha_star_const
y = (Z %*% alpha_star + d * beta_star) + errors[, 'epsilon']
# TEAM_IV Method
out = TEAM_IV(y, d, Z, K = 20, check_gam_sign = TRUE, verbose = 1L)
```
