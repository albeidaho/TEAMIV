# TEAMIV
Required functions for TEAM-IV method (Team-based Evaluation and Aggregation with MCP for Instrumental Variables)

Example usage:
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
