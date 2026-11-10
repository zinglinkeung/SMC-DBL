###############################################################################
#  realdata_analysis.R
###############################################################################

library(Rcpp); library(RcppArmadillo); library(glmnet)
library(quadprog); library(MASS); library(survival); library(smcfcs)
sourceCpp("iro_smc_fcs_lib.cpp")

# ================================================================
#  DATA PREPARATION
# ================================================================
df_raw <- read.csv("df_all.csv", stringsAsFactors = FALSE)

# Drop subjects with missing outcome
df_raw <- df_raw[!is.na(df_raw$OS_month) & !is.na(df_raw$DEAD), ]

# Identify SNP columns
snp_cols <- grep("^rs", colnames(df_raw), value = TRUE)

# Clinical predictors (all retained; missingness handled by imputation)
clinical_cols <- c("SEX", "AGE", "smoksort", "pc1", "pc2", "pc3",
                   "early_late", "RADS", "chemotx", "surgery")

predictor_cols <- c(snp_cols, clinical_cols)
p_loc <- length(predictor_cols)

# Build predictor matrix (numeric, with NAs preserved). SNPs are 0/1/2; the
# clinical indicators are read as numeric on their original 0/1 or 1/2 scales.
X_obs <- as.matrix(sapply(df_raw[, predictor_cols], as.numeric))
storage.mode(X_obs) <- "numeric"

# Standardise continuous columns. AGE and pc1-pc3 are continuous; SNPs and
# binary/ordinal clinical indicators are left on their original numeric scale,
# matching the convention in the simulation where SNP-style covariates were
# left on their additive scale.
cont_idx <- which(predictor_cols %in% c("AGE", "pc1", "pc2", "pc3"))
for(j in cont_idx) {
  mu <- mean(X_obs[, j], na.rm = TRUE)
  sd <- sd(X_obs[, j],  na.rm = TRUE)
  X_obs[, j] <- (X_obs[, j] - mu) / sd
}

# Outcome
Y     <- as.numeric(df_raw$OS_month)
delta <- as.numeric(df_raw$DEAD)

# Missingness indicator (1 = observed, 0 = missing)
R <- 1L * !is.na(X_obs)
storage.mode(R) <- "integer"

# Replace NAs in X_obs with NA (keep) for downstream functions
n_loc <- nrow(X_obs)

cat(sprintf("Data: n = %d, p = %d (%d SNPs + %d clinical)\n",
            n_loc, p_loc, length(snp_cols), length(clinical_cols)))
cat(sprintf("Events: %d (%.1f%% died)\n",
            sum(delta), 100 * mean(delta)))
cat(sprintf("Subjects with any missing predictor: %d (%.1f%%)\n",
            sum(rowSums(R == 0) > 0), 100 * mean(rowSums(R == 0) > 0)))
cat(sprintf("Cell-level missing rate: %.4f\n", mean(R == 0)))

# ================================================================
#  PARAMETERS
# ================================================================
L_outer_max <- 50
tol_outer   <- 1e-4
M_imp        <- 50  # paper Section 6 says M = 50 for the application
T_burn_chain <- 100 # paper Section 6 says T_0 = 100 for the application
S_inner      <- max(5, ceiling(log(n_loc * p_loc)))
nfold <- 10; n_lambda <- 100; tol_lasso <- 1e-6; maxiter_lasso <- 50000

cat(sprintf("Imputation params: M = %d, T_burn = %d, S_inner = %d\n",
            M_imp, T_burn_chain, S_inner))

# ================================================================
#  HELPERS 
# ================================================================
get_breslow <- function(Y, delta, X, beta)
  breslow_estimator_cpp(as.numeric(Y), as.numeric(delta), X, beta)

init_missing <- function(X_obs, R) {
  Xc <- X_obs
  for(j in 1:ncol(X_obs)) {
    mi <- which(R[, j] == 0); ob <- which(R[, j] == 1)
    if(length(mi) > 0 && length(ob) > 0)
      Xc[mi, j] <- sample(X_obs[ob, j], length(mi), replace = TRUE)
  }
  Xc
}

mean_impute <- function(X_obs, R) {
  Xc <- X_obs
  for(j in 1:ncol(X_obs)) {
    mi <- which(R[, j] == 0); ob <- which(R[, j] == 1)
    if(length(mi) > 0 && length(ob) > 0)
      Xc[mi, j] <- mean(X_obs[ob, j])
  }
  Xc
}

fit_cox_lasso <- function(X, Y, delta, lambda = NULL) {
  surv_y <- cbind(time = Y, status = delta)
  if(is.null(lambda)) {
    cvobj <- cv.glmnet(x = X, y = surv_y, family = "cox", alpha = 1,
                       standardize = FALSE, nfolds = nfold, nlambda = n_lambda)
    lambda <- cvobj$lambda.min
  }
  beta <- as.vector(coef(glmnet(x = X, y = surv_y, family = "cox", alpha = 1,
                                lambda = lambda, standardize = FALSE,
                                thresh = tol_lasso, maxit = maxiter_lasso)))
  list(beta = beta, lambda = lambda)
}

debias_at_multiplier <- function(X, Y, delta, beta_L1, neg_dll, ssq, r,
                                 multiplier) {
  n_l <- nrow(X); p_l <- ncol(X)
  gamma_n <- multiplier * sqrt(log(p_l) / n_l)
  my_pos <- which(r$values > 0); my_rank <- length(my_pos)
  beta_db <- rep(NA, p_l); se_db <- rep(NA, p_l)
  Theta <- matrix(NA, p_l, p_l)
  if(my_rank > 0) {
    Dmat <- diag(r$values[my_pos], nrow = my_rank)
    dvec <- rep(0, my_rank)
    Amat <- t(rbind(-r$vectors[, my_pos] %*% Dmat,
                     r$vectors[, my_pos] %*% Dmat))
    for(j in 1:p_l) {
      e_j <- rep(0, p_l); e_j[j] <- 1
      bvec <- c(-e_j, e_j) - gamma_n * rep(1, 2 * p_l)
      res <- tryCatch(solve.QP(Dmat = Dmat, dvec = dvec,
                               Amat = Amat, bvec = bvec),
                      error = function(e) NULL)
      if(is.null(res)) {
        m <- solve(ssq + gamma_n * diag(p_l), e_j)
      } else {
        m <- as.vector(r$vectors[, my_pos] %*% res$solution)
        if(my_rank < p_l)
          m <- m + as.vector(r$vectors[, -my_pos, drop = FALSE] %*%
                             rep(0, p_l - my_rank))
      }
      beta_db[j] <- beta_L1[j] - sum(m * neg_dll)
      mSm <- as.numeric(t(m) %*% ssq %*% m)
      se_db[j] <- sqrt(max(mSm, 1e-12) / n_l)
      Theta[j, ] <- m
    }
  } else {
    Theta <- solve(ssq + 0.1 * diag(p_l))
    beta_db <- as.vector(beta_L1 - Theta %*% neg_dll)
    var_diag <- diag(Theta %*% ssq %*% t(Theta)) / n_l
    se_db <- sqrt(pmax(var_diag, 1e-12))
  }
  list(beta_db = beta_db, se = se_db, Theta = Theta)
}

debiased_lasso_cox <- function(X, Y, delta, multiplier = 0.5) {
  n_l <- nrow(X); p_l <- ncol(X)
  surv_y <- cbind(time = Y, status = delta)
  cvobj <- cv.glmnet(x = X, y = surv_y, family = "cox", alpha = 1,
                     standardize = FALSE, nfolds = nfold, nlambda = n_lambda)
  beta_L1 <- as.vector(coef(glmnet(x = X, y = surv_y, family = "cox", alpha = 1,
                                   lambda = cvobj$lambda.min,
                                   standardize = FALSE, thresh = tol_lasso,
                                   maxit = maxiter_lasso)))
  neg_ll <- 0; neg_dll <- rep(0, p_l)
  neg_ddll <- matrix(0, p_l, p_l); ssq <- matrix(0, p_l, p_l)
  neg_loglik_functions_cpp_ext(neg_ll, neg_dll, neg_ddll, ssq,
                               X, Y, delta, beta_L1)
  r <- eigen(ssq); r$values[r$values <= 1e-14] <- 0
  db <- debias_at_multiplier(X, Y, delta, beta_L1, neg_dll, ssq, r, multiplier)
  list(beta_db = db$beta_db, se = db$se, Theta = db$Theta,
       beta_L1 = beta_L1, ssq = ssq, lambda = cvobj$lambda.min)
}

tune_ridge_multiplier <- function(Xc, R, miss_vars, n_loc, p_loc, n_cv = 5) {
  j <- miss_vars[1] + 1L
  obs_idx <- which(R[, j] == 1)
  n_obs <- length(obs_idx); pm1 <- p_loc - 1
  y_obs <- Xc[obs_idx, j]
  Z_obs <- Xc[obs_idx, -j, drop = FALSE]
  base_rate <- p_loc * log(max(p_loc, 2)) / n_loc
  b_upper <- max(2, 20 * sqrt(p_loc / n_loc))
  b_candidates <- exp(seq(from = log(0.01), to = log(b_upper),
                          length.out = 15))
  cv_idx <- sample(rep(1:n_cv, length.out = n_obs))
  all_cv_mse <- rep(0, length(b_candidates))
  for(bb in seq_along(b_candidates)) {
    b_val <- b_candidates[bb]; lam_c <- b_val * base_rate
    mse_per_fold <- numeric(n_cv)
    for(k in 1:n_cv) {
      tr <- which(cv_idx != k); te <- which(cv_idx == k)
      yt <- y_obs[tr]; ye <- y_obs[te]
      Zt <- Z_obs[tr, , drop = FALSE]; Ze <- Z_obs[te, , drop = FALSE]
      ym <- mean(yt); zm <- colMeans(Zt)
      yc <- yt - ym; Zc <- sweep(Zt, 2, zm); Zce <- sweep(Ze, 2, zm)
      ZtZ <- crossprod(Zc) / length(tr); diag(ZtZ) <- diag(ZtZ) + lam_c
      Zty <- crossprod(Zc, yc) / length(tr)
      gh <- tryCatch(solve(ZtZ, Zty), error = function(e) rep(0, pm1))
      mse_per_fold[k] <- mean((ye - (ym + as.vector(Zce %*% gh)))^2)
    }
    all_cv_mse[bb] <- mean(mse_per_fold)
  }
  best_b <- b_candidates[which.min(all_cv_mse)]
  cat(sprintf("    Tuned ridge multiplier: b = %.4f\n", best_b))
  best_b
}

permute_sweep <- function(miss_vars, cv_lam) {
  perm <- sample.int(length(miss_vars))
  list(mv = miss_vars[perm], cl = cv_lam[perm])
}

# ================================================================
#  METHOD 1: COMPLETE-CASE DEBIASED LASSO  
# ================================================================
complete_case_dl <- function(X_obs, Y, delta, R) {
  cc <- which(rowSums(R == 0) == 0)
  cat(sprintf("    Complete-case rows: %d / %d\n", length(cc), nrow(X_obs)))
  Xcc <- X_obs[cc, , drop = FALSE]; Ycc <- Y[cc]; dcc <- delta[cc]
  r <- debiased_lasso_cox(Xcc, Ycc, dcc)
  list(beta_hat = r$beta_db, se = r$se,
       ci_lo = r$beta_db - 1.96 * r$se,
       ci_hi = r$beta_db + 1.96 * r$se,
       n_used = length(cc))
}

# ================================================================
#  METHOD 2: SMC-DBL  
# ================================================================
smc_dbl <- function(X_obs, Y, delta, R) {
  n_loc <- nrow(X_obs); p_loc <- ncol(X_obs)
  miss_vars <- as.integer(which(colSums(R == 0) > 0) - 1L)
  if(length(miss_vars) == 0) stop("No missing variables")
  lambda_j <- 0.5 * sqrt(log(max(p_loc, 2)) / n_loc)
  R_int <- matrix(as.integer(R), n_loc, p_loc)
  n_mv <- length(miss_vars)

  # Phase 0: tuning chain
  Xc <- init_missing(X_obs, R)
  cx <- fit_cox_lasso(Xc, Y, delta)
  beta_cur <- cx$beta; lambda_n <- cx$lambda
  bres <- get_breslow(Y, delta, Xc, beta_cur)

  base_rate <- p_loc * log(max(p_loc, 2)) / n_loc
  cv_lam_init <- rep(0.1 * base_rate, n_mv)

  cat("    Phase 0: tuning chain ...\n")
  for(l in 1:L_outer_max) {
    for(s in 1:S_inner) {
      ps <- permute_sweep(miss_vars, cv_lam_init)
      smc_fcs_sweep_ridge_bayes_cpp(Xc, Y, delta, R_int, beta_cur,
                                    bres$times, bres$values,
                                    lambda_j, ps$mv, ps$cl)
    }
    beta_new <- as.vector(coef(glmnet(
      x = Xc, y = cbind(time = Y, status = delta), family = "cox", alpha = 1,
      lambda = lambda_n, standardize = FALSE,
      thresh = tol_lasso, maxit = maxiter_lasso)))
    bres <- get_breslow(Y, delta, Xc, beta_new)
    if(sum(abs(beta_new - beta_cur)) < tol_outer) {
      cat(sprintf("    Phase 0 converged at iter %d\n", l))
      beta_cur <- beta_new; break
    }
    beta_cur <- beta_new
  }

  qp_multiplier <- 0.5
  ridge_b <- tune_ridge_multiplier(Xc, R, miss_vars, n_loc, p_loc)
  cv_lam <- rep(ridge_b * base_rate, n_mv)

  # Phase 1+2: M independent inferential chains
  run_one_chain <- function(seed) {
    set.seed(seed)
    Xc_m <- init_missing(X_obs, R)
    cx_m <- fit_cox_lasso(Xc_m, Y, delta)
    beta_m <- cx_m$beta; lambda_n_m <- cx_m$lambda
    bres_m <- get_breslow(Y, delta, Xc_m, beta_m)
    for(l in 1:L_outer_max) {
      for(s in 1:S_inner) {
        ps <- permute_sweep(miss_vars, cv_lam_init)
        smc_fcs_sweep_ridge_bayes_cpp(Xc_m, Y, delta, R_int, beta_m,
                                      bres_m$times, bres_m$values,
                                      lambda_j, ps$mv, ps$cl)
      }
      beta_new_m <- as.vector(coef(glmnet(
        x = Xc_m, y = cbind(time = Y, status = delta), family = "cox", alpha = 1,
        lambda = lambda_n_m, standardize = FALSE,
        thresh = tol_lasso, maxit = maxiter_lasso)))
      bres_m <- get_breslow(Y, delta, Xc_m, beta_new_m)
      if(sum(abs(beta_new_m - beta_m)) < tol_outer) {
        beta_m <- beta_new_m; break
      }
      beta_m <- beta_new_m
    }
    ridge_cache <- precompute_ridge_cache_cpp(Xc_m, R_int, miss_vars, cv_lam)
    for(t in 1:T_burn_chain) {
      for(s in 1:S_inner) {
        ps <- permute_sweep(miss_vars, cv_lam)
        smc_fcs_sweep_ridge_bayes_cached_cpp(
          Xc_m, Y, delta, R_int, beta_m,
          bres_m$times, bres_m$values,
          lambda_j, ps$mv, ridge_cache)
      }
      bres_m <- get_breslow(Y, delta, Xc_m, beta_m)
    }
    list(Xc = Xc_m)
  }

  cat(sprintf("    Phase 1+2: %d chains ...\n", M_imp))
  set.seed(20251104)
  chain_seeds <- sample.int(1e7, M_imp)
  chain_results <- vector("list", M_imp)
  for(m in 1:M_imp) {
    cat(sprintf("      chain %d/%d\n", m, M_imp))
    chain_results[[m]] <- run_one_chain(chain_seeds[m])
  }

  # Phase 3: per-chain debiased lasso + Rubin pooling
  cat("    Phase 3: per-chain debiased lasso + Rubin pooling ...\n")
  bdb <- matrix(NA, M_imp, p_loc)
  W_per_chain <- matrix(NA, M_imp, p_loc)
  for(m in 1:M_imp) {
    Xm <- chain_results[[m]]$Xc
    r <- debiased_lasso_cox(Xm, Y, delta, multiplier = qp_multiplier)
    bdb[m, ] <- r$beta_db; W_per_chain[m, ] <- r$se^2
  }
  bb <- colMeans(bdb)
  VW <- colMeans(W_per_chain)
  VB <- apply(bdb, 2, var)
  V_total <- VW + (1 + 1/M_imp) * VB
  se <- sqrt(pmax(V_total, 1e-12))
  rm_ratio <- (1 + 1/M_imp) * VB / pmax(VW, 1e-12)
  df_rubin <- pmax((M_imp - 1) * (1 + 1/rm_ratio)^2, 2)
  list(beta_hat = bb, se = se,
       ci_lo = bb - qt(0.975, df_rubin) * se,
       ci_hi = bb + qt(0.975, df_rubin) * se,
       VB = VB, VW = VW)
}

# ================================================================
#  METHOD 3: STANDARD IRO 
# ================================================================
standard_iro <- function(X_obs, Y, delta, R) {
  n_loc <- nrow(X_obs); p_loc <- ncol(X_obs)
  miss_vars <- as.integer(which(colSums(R == 0) > 0) - 1L)
  lambda_j <- 0.5 * sqrt(log(max(p_loc, 2)) / n_loc)
  R_int <- matrix(as.integer(R), n_loc, p_loc)
  T_total <- max(20, ceiling(3 * log(n_loc * p_loc))) +
             M_imp * max(10, ceiling(2 * log(n_loc * p_loc)))
  T_burn  <- max(20, ceiling(3 * log(n_loc * p_loc)))

  Xc <- init_missing(X_obs, R)
  cx <- fit_cox_lasso(Xc, Y, delta)
  beta_cur <- cx$beta; lambda_n <- cx$lambda
  bres <- get_breslow(Y, delta, Xc, beta_cur)

  beta_chain <- matrix(NA, T_total, p_loc)
  for(t in 1:T_total) {
    for(s in 1:S_inner)
      smc_fcs_sweep_cpp(Xc, Y, delta, R_int, beta_cur,
                        bres$times, bres$values, lambda_j, miss_vars)
    if(t %% 10 == 1) {
      cx <- fit_cox_lasso(Xc, Y, delta); lambda_n <- cx$lambda
    }
    beta_cur <- as.vector(coef(glmnet(
      x = Xc, y = cbind(time = Y, status = delta), family = "cox", alpha = 1,
      lambda = lambda_n, standardize = FALSE,
      thresh = tol_lasso, maxit = maxiter_lasso)))
    bres <- get_breslow(Y, delta, Xc, beta_cur)
    if(t > T_burn) beta_chain[t - T_burn, ] <- beta_cur
  }
  n_coll <- T_total - T_burn
  beta_chain <- beta_chain[1:n_coll, , drop = FALSE]
  beta_bar <- colMeans(beta_chain)

  neg_ll <- 0; neg_dll <- rep(0, p_loc)
  neg_ddll <- matrix(0, p_loc, p_loc); ssq <- matrix(0, p_loc, p_loc)
  neg_loglik_functions_cpp_ext(neg_ll, neg_dll, neg_ddll, ssq,
                               Xc, Y, delta, beta_cur)
  hess_inv <- tryCatch(solve(neg_ddll + 0.01 * diag(p_loc)),
                       error = function(e) solve(ssq + 0.1 * diag(p_loc)))
  VW <- pmax(diag(hess_inv), 1e-12) / n_loc
  VB <- apply(beta_chain, 2, var)
  V_total <- VW + VB / n_coll
  se <- sqrt(pmax(V_total, 1e-12))

  list(beta_hat = beta_bar, se = se,
       ci_lo = beta_bar - 1.96 * se,
       ci_hi = beta_bar + 1.96 * se)
}

# ================================================================
#  METHOD 4: STANDARD SMC-FCS + DL + RUBIN
# ================================================================
standard_smc_fcs <- function(X_obs, Y, delta, R) {
  n_loc <- nrow(X_obs); p_loc <- ncol(X_obs)
  df <- as.data.frame(X_obs)
  colnames(df) <- paste0("X", 1:p_loc)
  df$time <- Y; df$status <- as.integer(delta)
  xvars <- paste0("X", 1:p_loc)
  sm_formula <- paste0("Surv(time, status) ~ ", paste(xvars, collapse = " + "))
  method_vec <- rep("", ncol(df)); names(method_vec) <- colnames(df)
  miss_cols <- which(colSums(R == 0) > 0)
  for(j in miss_cols) method_vec[paste0("X", j)] <- "norm"
  imps <- smcfcs::smcfcs(originaldata = df, smtype = "coxph",
                         smformula = sm_formula, method = method_vec,
                         m = M_imp, numit = 20, rjlimit = 1000, noisy = FALSE)
  bdb <- matrix(NA, M_imp, p_loc); Wpc <- matrix(NA, M_imp, p_loc)
  for(m in 1:M_imp) {
    Xm <- as.matrix(imps$impDatasets[[m]][, xvars])
    r <- debiased_lasso_cox(Xm, Y, delta)
    bdb[m, ] <- r$beta_db; Wpc[m, ] <- r$se^2
  }
  bb <- colMeans(bdb); VW <- colMeans(Wpc); VB <- apply(bdb, 2, var)
  V_total <- VW + (1 + 1/M_imp) * VB
  se <- sqrt(pmax(V_total, 1e-12))
  rm_ratio <- (1 + 1/M_imp) * VB / pmax(VW, 1e-12)
  df_rubin <- pmax((M_imp - 1) * (1 + 1/rm_ratio)^2, 2)
  list(beta_hat = bb, se = se,
       ci_lo = bb - qt(0.975, df_rubin) * se,
       ci_hi = bb + qt(0.975, df_rubin) * se)
}

# ================================================================
#  METHOD 5: MEAN IMPUTATION + DEBIASED LASSO
# ================================================================
mean_imp_dl <- function(X_obs, Y, delta, R) {
  Xc <- mean_impute(X_obs, R)
  r <- debiased_lasso_cox(Xc, Y, delta)
  list(beta_hat = r$beta_db, se = r$se,
       ci_lo = r$beta_db - 1.96 * r$se,
       ci_hi = r$beta_db + 1.96 * r$se)
}

# ================================================================
#  RUN ALL METHODS
# ================================================================
methods <- c("complete_case", "smc_dbl", "standard_iro",
             "std_smc_fcs", "mean_imp")
labels  <- c("Complete-case DL", "SMC-DBL", "Standard IRO (Liang)",
             "Std SMC-FCS", "Mean Imp + DL")
method_fns <- list(
  complete_case = function() complete_case_dl(X_obs, Y, delta, R),
  smc_dbl       = function() smc_dbl(X_obs, Y, delta, R),
  standard_iro  = function() standard_iro(X_obs, Y, delta, R),
  std_smc_fcs   = function() standard_smc_fcs(X_obs, Y, delta, R),
  mean_imp      = function() mean_imp_dl(X_obs, Y, delta, R)
)

res_all <- list()
timings <- list()
for(method in methods) {
  cat(sprintf("\n>>> Running %s ...\n", method))
  t0 <- proc.time()
  res <- tryCatch(method_fns[[method]](),
                  error = function(e) {
                    cat(sprintf("    FAILED: %s\n", e$message)); NULL
                  })
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("    Done (%.1fs)\n", elapsed))
  if(!is.null(res)) res_all[[method]] <- res
  timings[[method]] <- elapsed
}

# ================================================================
#  SAVE
# ================================================================
dir.create("realdata_results", showWarnings = FALSE)
save(res_all, timings, predictor_cols, snp_cols, clinical_cols,
     n_loc, p_loc, file = "realdata_results/results.RData")

# Tidy summary CSV: one row per (method, predictor)
tidy <- data.frame()
for(mi in seq_along(methods)) {
  mn <- methods[mi]; ml <- labels[mi]
  if(is.null(res_all[[mn]])) next
  r <- res_all[[mn]]
  is_snp <- predictor_cols %in% snp_cols
  z <- r$beta_hat / pmax(r$se, 1e-12)
  pval <- 2 * pnorm(-abs(z))
  qval_snp <- rep(NA, length(predictor_cols))
  qval_snp[is_snp] <- p.adjust(pval[is_snp], method = "BH")
  tidy <- rbind(tidy, data.frame(
    method   = ml,
    variable = predictor_cols,
    is_snp   = is_snp,
    beta_hat = r$beta_hat,
    se       = r$se,
    ci_lo    = r$ci_lo,
    ci_hi    = r$ci_hi,
    z        = z,
    pval     = pval,
    qval_snp_BH = qval_snp
  ))
}
write.csv(tidy, "realdata_results/summary.csv", row.names = FALSE)

cat("\nWrote realdata_results/results.RData and realdata_results/summary.csv\n")
cat(sprintf("Total runtime: %.1f minutes\n",
            sum(unlist(timings)) / 60))
