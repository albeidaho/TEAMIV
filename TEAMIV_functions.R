library(dplyr)
library(tibble)
library(purrr)
library(glmnet)
library(ncvreg)

cv_loss_teamiv <- function(y, d, Z, indices_invalid = integer(0),
                           K = 10, seed = NULL, folds = NULL) {
  # Purpose: compute K-fold CV errors for the TEAM-IV-style linear model
  #          y ~ 1 + P_Z d + Z_{A0},
  #          i.e., least-squares on y with regressors (1, P_Z d, Z_invalid).
  #
  # Returns:
  #   list(
  #     mean_mse, se_mse,
  #     mean_rmse, se_rmse,
  #     fold_mse, fold_rmse,
  #     folds
  #   )

  n <- length(y)
  if (length(d) != n || nrow(Z) != n) {
    stop("Lengths of y, d and number of rows of Z must all match.")
  }

  ## Folds
  if (is.null(folds)) {
    if (!is.null(seed)) set.seed(seed)
    folds <- sample(rep(1:K, length.out = n))
  } else {
    if (length(folds) != n) stop("folds must have length equal to n")
    if (!all(folds %in% 1:K)) stop("folds must contain only integers 1..K")
  }

  fold_mse  <- numeric(K)
  fold_rmse <- numeric(K)

  for (k in 1:K) {
    test_idx  <- which(folds == k)
    train_idx <- setdiff(seq_len(n), test_idx)

    y_train <- y[train_idx]
    d_train <- d[train_idx]
    Z_train <- Z[train_idx, , drop = FALSE]

    y_test <- y[test_idx]
    d_test <- d[test_idx]
    Z_test <- Z[test_idx, , drop = FALSE]

    ## Projection matrices per fold
    # P_Z = Z (Z^T Z)^{-1} Z^T, implemented via qr.solve for stability
    PZ_train <- Z_train %*% qr.solve(t(Z_train) %*% Z_train) %*% t(Z_train)
    PZ_test  <- Z_test  %*% qr.solve(t(Z_test)  %*% Z_test)  %*% t(Z_test)

    ## Projected d for TEAM-IV design
    d_proj_train <- PZ_train %*% d_train
    d_proj_test  <- PZ_test  %*% d_test

    ## Invalid instruments (can be empty)
    Zinv_train <- if (length(indices_invalid) > 0L) {
      Z_train[, indices_invalid, drop = FALSE]
    } else {
      NULL
    }

    Zinv_test <- if (length(indices_invalid) > 0L) {
      Z_test[, indices_invalid, drop = FALSE]
    } else {
      NULL
    }

    ## Design matrices
    X_train <- cbind(1, d_proj_train, Zinv_train)
    X_test  <- cbind(1, d_proj_test,  Zinv_test)

    ## Fit TEAM-IV-style linear model (unpenalized)
    fit      <- stats::lm.fit(x = X_train, y = y_train)
    coef_hat <- fit$coefficients

    ## Predicted y on test set
    yhat_test <- c(X_test %*% coef_hat)

    ## Fold-level prediction error for y (TEAM-IV loss)
    resid_test <- y_test - yhat_test
    fold_mse[k]  <- mean(resid_test^2)
    fold_rmse[k] <- sqrt(fold_mse[k])
  }

  ## Summaries over folds
  mean_mse     <- mean(fold_mse)
  se_mean_mse  <- sqrt(var(fold_mse) / K)
  mean_rmse    <- sqrt(mean_mse)
  se_rmse      <- se_mean_mse / (2 * mean_rmse)

  list(
    mean_mse  = mean_mse,
    se_mse    = se_mean_mse,
    mean_rmse = mean_rmse,
    se_rmse   = se_rmse,
    fold_mse  = fold_mse,
    fold_rmse = fold_rmse,
    folds     = folds
  )
}

mcp_valid_team_sign2 <- function(
  y,
  d,
  Z,
  lam_se = c(1, 0.5, 0.25, 0.1),
  max_valid_size = Inf,
  alpha_tol = 1e-6,
  K = 10
) {
  stopifnot(is.matrix(Z), is.numeric(y), is.numeric(d))
  n <- length(y)
  p <- ncol(Z)

  y <- as.numeric(y)
  d <- as.numeric(d)
  Z <- as.matrix(Z)

  # increment amount per lam_se value
  inc_amt <- 1 / length(lam_se)

  # initialize symmetric "teamness" matrix
  results <- matrix(0, nrow = p, ncol = p)
  rownames(results) <- colnames(results) <- paste0("Z", seq_len(p))

  folds0 <- ncvreg::assign_fold(y, K)

  # helper: CV 1-SE rule for arbitrary se
  lambda_xse <- function(cvfit, se) {
    cvtib <- with(cvfit, tibble::tibble(cve = cve, cvse = cvse, lam = lambda))
    if (nrow(cvtib) == 0) {
      return(NA_real_)
    }
    thr <- min(cvtib$cve) + se * cvtib$cvse[which.min(cvtib$cve)]
    sub <- dplyr::filter(cvtib, cve <= thr)
    if (nrow(sub) == 0) {
      return(NA_real_)
    }
    dplyr::arrange(sub, dplyr::desc(lam)) %>%
      dplyr::slice(1) %>%
      dplyr::pull(lam)
  }

  # loop over odd contiguous IV set sizes
  for (valid_size in seq(3, p, by = 2L)) {
    if (valid_size > max_valid_size) {
      break
    }
    nsets <- p - valid_size + 1L

    for (start in seq_len(nsets)) {
      considered_inds <- seq.int(start, start + valid_size - 1L)

      # projection PZ
      PZ <- Z[, considered_inds, drop = FALSE] %*%
        qr.solve(t(Z[, considered_inds]) %*% Z[, considered_inds]) %*%
        t(Z[, considered_inds])

      d_proj <- PZ %*% d
      X_proj <- cbind(d_proj, Z[, considered_inds, drop = FALSE])
      penalty_factor <- c(0, rep(1, ncol(X_proj) - 1))

      # fit sisVIVE-style MCP model
      cvncvreg_out <- tryCatch(
        ncvreg::cv.ncvreg(
          X_proj,
          y,
          penalty = "MCP",
          penalty.factor = penalty_factor,
          fold = folds0
        ),
        error = function(e) NULL
      )
      if (is.null(cvncvreg_out)) {
        next
      }

      # ---- LOOP OVER lam_se VECTOR ----
      for (se_val in lam_se) {
        lam0 <- lambda_xse(cvncvreg_out, se = se_val)
        cvncvreg_alpha_hat <- coef(cvncvreg_out, lam = lam0)[-(1:2)]
        if (is.null(cvncvreg_alpha_hat)) {
          next
        }

        sisvive_alpha_results_binary <- rep(NA, p)
        alphas <- cvncvreg_alpha_hat

        # valid = 0, invalid = ±1
        sisvive_alpha_results_binary[considered_inds] <- ifelse(
          abs(alphas) <= alpha_tol,
          0,
          sign(alphas)
        )

        valid_inds <- which(sisvive_alpha_results_binary == 0)
        neginvalid_inds <- which(sisvive_alpha_results_binary == -1)
        posinvalid_inds <- which(sisvive_alpha_results_binary == 1)
        invalid_inds <- union(neginvalid_inds, posinvalid_inds)

        # increment: both valid
        if (length(valid_inds) > 1) {
          results[valid_inds, valid_inds] <-
            results[valid_inds, valid_inds] + inc_amt
        }

        # decrement: mixed valid/invalid
        if (length(valid_inds) > 0 && length(invalid_inds) > 0) {
          results[valid_inds, invalid_inds] <-
            results[valid_inds, invalid_inds] - inc_amt
          results[invalid_inds, valid_inds] <-
            results[invalid_inds, valid_inds] - inc_amt
        }

        if (length(neginvalid_inds) > 0 && length(posinvalid_inds) > 0) {
          # decrement negative–positive invalid pairings
          results[neginvalid_inds, posinvalid_inds] <-
            results[neginvalid_inds, posinvalid_inds] - inc_amt
          results[posinvalid_inds, neginvalid_inds] <-
            results[posinvalid_inds, neginvalid_inds] - inc_amt
        }
      } # end lam_se loop
    } # end start loop
  } # end valid_size loop

  # summarize into tibble
  results_tib <- as_tibble(results, .name_repair = "minimal") %>%
    mutate(iv = paste0("Z", seq_len(p))) %>%
    rowwise() %>%
    mutate(
      same_team_inds = list(which(c_across(starts_with("Z")) > 0)),
      same_team = list(paste0("Z", same_team_inds))
    ) %>%
    ungroup() %>%
    dplyr::select(iv, same_team_inds)

  results_tib_unique <- results_tib %>%
    dplyr::filter(lengths(same_team_inds) > 1) %>%
    distinct(same_team_inds, .keep_all = TRUE) %>%
    dplyr::arrange(desc(lengths(same_team_inds)))

  return(list(res_tib_unique = results_tib_unique, teamness_mat = results))
}

TEAM_IV <- function(
  y,
  d,
  Z,
  K = 20,
  lam_se = 0.1,
  check_gam_sign = FALSE,
  verbose = 0L,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  folds0 <- ncvreg::assign_fold(y, K, seed)

  nZ <- ncol(Z)
  if (nZ < 3) {
    stop("Need at least 3 instruments")
  }
# ---- Step 0. Adjust exposure for average gamma_hat
  ridge_out_pre <- glmnet::cv.glmnet(Z, d, alpha = 0)
  gamma_hat <- as.numeric(coef(ridge_out_pre, s = 'lambda.min'))[-1]
  if(check_gam_sign && any(gamma_hat<0)){
    Z <- sweep(Z, MARGIN = 2, STATS = as.integer(gamma_hat<0)*(-2)+1, FUN = "*")
    ridge_out_post <- glmnet::cv.glmnet(Z, d, alpha = 0)
    gamma_hat <- as.numeric(coef(ridge_out_post, s = 'lambda.min'))[-1]
  }
  mean_gamma_hat <- mean(gamma_hat)
  d_adj <- as.numeric(d - Z %*% (gamma_hat - mean_gamma_hat))

  # ---- Step 1. Ordering from check_equals
  PZ = Z %*%
    qr.solve(
      t(Z) %*%
        Z
    ) %*%
    t(Z)
  d_proj = PZ %*% d
  X <- cbind(d_proj, Z)
  penalty_factor <- c(0, rep(1, ncol(X) - 1))
  res_ridge <- glmnet::cv.glmnet(
    X,
    y,
    penalty.factor = penalty_factor,
    alpha = 0
  )
  alpha_hat = coef(res_ridge, s = res_ridge$lambda.min) %>%
    as.matrix %>%
    .[-(1:2), ]
  names(alpha_hat) = paste0(1:length(alpha_hat))
  ord <- alpha_hat %>%
    as.data.frame %>%
    tibble::rownames_to_column() %>%
    dplyr::rename(alpha_hat = 2) %>%
    dplyr::arrange(alpha_hat) %>%
    dplyr::pull(rowname) %>%
    as.integer
  inv_ord = sapply(1:nZ, function(x) which(ord==x))
  Z_ord <- Z[, ord, drop = FALSE]

  # helper: obtain sparsest cv.ncvreg() lambda within certain SE threshold
  lambda_xse <- function(cvfit, se = 1) {
    cvtib <- with(cvfit, tibble::tibble(cve = cve, cvse = cvse, lam = lambda))
    if (nrow(cvtib) == 0) {
      return(NA_real_)
    }
    thr <- min(cvtib$cve) + se * cvtib$cvse[which.min(cvtib$cve)]
    sub <- cvtib %>% dplyr::filter(cve <= thr)
    if (nrow(sub) == 0) {
      return(NA_real_)
    }
    sub %>%
      dplyr::arrange(dplyr::desc(lam)) %>%
      dplyr::slice(1) %>%
      dplyr::pull(lam)
  }

  # ---- Helper: evaluate a candidate valid IV set
  evaluate_candidate <- function(valid_inds, folds = folds0) {
    if (length(valid_inds) < 2) {
      return(list(
        mean_rmse = Inf,
        se_rmse = Inf,
        fold_mse = rep(Inf, K),
        beta_hat = NA
      ))
    }
    invalid_inds <- setdiff(seq_len(nZ), valid_inds)
    cvsl <- tryCatch(
      cv_loss_teamiv(
        y,
        d,
        Z_ord,
        indices_invalid = invalid_inds,
        K = K,
        folds = folds
      ),
      error = function(e) NULL
    )

    beta_hat <- tryCatch(
      {
        PZ = Z_ord[, valid_inds, drop = FALSE] %*%
          qr.solve(
            t(Z_ord[, valid_inds, drop = FALSE]) %*%
              Z_ord[, valid_inds, drop = FALSE]
          ) %*%
          t(Z_ord[, valid_inds, drop = FALSE])
        d_proj = PZ %*% d
        X <- cbind(d_proj, Z_ord[, valid_inds, drop = FALSE])
        penalty_factor <- c(0, rep(1, ncol(X) - 1))
        ##
        res_mcp <- ncvreg::cv.ncvreg(
          X,
          y,
          penalty = "MCP",
          penalty.factor = penalty_factor,
          alpha = 1
        )
        lam0_mcp = lambda_xse(res_mcp, se = lam_se)
        ##
        res_50 <- ncvreg::cv.ncvreg(
          X,
          y,
          penalty.factor = penalty_factor,
          alpha = 0.50
        )
        min_cve_a50 = min(res_50$cve)
        lam0_50 = lambda_xse(res_50, se = lam_se)
        ##
        res_75 <- ncvreg::cv.ncvreg(
          X,
          y,
          penalty.factor = penalty_factor,
          alpha = 0.75
        )
        min_cve_a75 = min(res_75$cve)
        lam0_75 = lambda_xse(res_75, se = lam_se)
        min_a_val = which.min(c(min_cve_a50, min_cve_a75))
        if (!is.null(res_50) & !is.null(lam0_50)) {
          bhat_res_50 = coef(res_50, lam = lam0_50)[2]
          bhat_res_75 = coef(res_75, lam = lam0_75)[2]
          bhat_min_mn = c(bhat_res_50, bhat_res_75)[min_a_val]
          c(
            bhat_min_mn,
            bhat_res_50,
            bhat_res_75,
            coef(res_mcp, lam = lam0_mcp)[2]
          ) %>%
            return
        }
      },
      error = function(e) NA
    )

    mean_rmse <- if (!is.null(cvsl)) cvsl$mean_rmse else Inf
    se_rmse <- if (!is.null(cvsl)) cvsl$se_rmse else Inf
    fold_mse <- if (!is.null(cvsl)) cvsl$fold_mse else rep(Inf, K)

    list(
      mean_rmse = mean_rmse,
      se_rmse = se_rmse,
      fold_mse = fold_mse,
      beta_afl = beta_hat[1],
      beta_a50 = beta_hat[2],
      beta_a75 = beta_hat[3],
      beta_mcp = beta_hat[4]
    )
  }

  # ---- Helper: union of teams within x SEs of best
  union_xse_valid_set <- function(pair_eval, se_num = 2) {
    if (nrow(pair_eval) == 0) {
      return(NULL)
    }
    min_rmse_row <- pair_eval %>% dplyr::arrange(mean_rmse) %>% slice(1)
    threshold <- min_rmse_row$mean_rmse + se_num * min_rmse_row$se_rmse

    cluster_rows <- pair_eval %>% dplyr::filter(mean_rmse <= threshold)
    flat_positions <- sort(unique(unlist(cluster_rows$valid_ivs)))
    union_valid_inds <- inv_ord[flat_positions]
    eval_union <- evaluate_candidate(union_valid_inds)

    tibble(
      iv = paste0('union_', se_num, 'se'),
      valid_ivs = list(flat_positions),
      n_valid = length(flat_positions),
      mean_rmse = eval_union$mean_rmse,
      se_rmse = eval_union$se_rmse,
      beta_afl = as.numeric(eval_union$beta_afl),
      beta_a50 = as.numeric(eval_union$beta_a50),
      beta_a75 = as.numeric(eval_union$beta_a75),
      beta_mcp = as.numeric(eval_union$beta_mcp),
      fold_mse = list(eval_union$fold_mse)
    )
  }

  # ---- Step 3. Identify sisVIVE valid teams
  mcp_valid_team_out = mcp_valid_team_sign2(
    y,
    d_adj,
    Z_ord,
    lam_se = c(1, 0.5, 0.25, 0.1),
    K=K
  )

  # ---- Step 4. Evaluate each team
  pair_eval <- mcp_valid_team_out$res_tib_unique %>%
    rowwise() %>%
    mutate(
      eval = list(evaluate_candidate(same_team_inds)),
      n_valid = length(same_team_inds),
      fold_mse = list(eval$fold_mse),
      original_iv_ind = list(ord[unlist(same_team_inds)])
    ) %>%
    ungroup() %>%
    mutate(
      mean_rmse = purrr::map_dbl(eval, ~ .x$mean_rmse),
      se_rmse = purrr::map_dbl(eval, ~ .x$se_rmse),
      beta_afl = purrr::map_dbl(eval, ~ as.numeric(.x$beta_afl)),
      beta_a50 = purrr::map_dbl(eval, ~ as.numeric(.x$beta_a50)),
      beta_a75 = purrr::map_dbl(eval, ~ as.numeric(.x$beta_a75)),
      beta_mcp = purrr::map_dbl(eval, ~ as.numeric(.x$beta_mcp))
    ) %>%
    dplyr::select(
      iv,
      valid_ivs = original_iv_ind,
      n_valid,
      mean_rmse,
      se_rmse,
      beta_afl,
      beta_a50,
      beta_a75,
      beta_mcp,
      fold_mse
    )

  # ---- Step 5. Compute unions within 1, 2, 5 SE of min
  unions <- list(1, 2, 5) %>%
    purrr::map(~ union_xse_valid_set(pair_eval, se_num = .x)) %>%
    purrr::compact() %>%
    dplyr::bind_rows()

  # ---- Step 6. Final combined output
  return_tib <- dplyr::bind_rows(pair_eval, unions)

  if(verbose>=1L){print(mcp_valid_team_out$teamness_mat); print(return_tib)}

  return(return_tib)
}