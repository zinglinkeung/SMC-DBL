###############################################################################
#  bench_single_iter.R  —  M independent serial chains
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])


library(Rcpp); library(RcppArmadillo); library(glmnet)
library(quadprog); library(MASS); library(survival); library(smcfcs)
sourceCpp("iro_smc_fcs_lib.cpp")

# ================================================================
#  SETTINGS TABLE
# ================================================================
settings <- expand.grid(
  n   = c(500, 1000, 2000),
  p   = c(20, 50, 100, 200)
)
settings$s0 <- 5
n_settings <- nrow(settings)
n_reps     <- 100

setting_id <- ((task_id - 1) %% n_settings) + 1
rep_id     <- ((task_id - 1) %/% n_settings) + 1

n  <- settings$n[setting_id]
p  <- settings$p[setting_id]
s0 <- settings$s0[setting_id]
rho <- 0.5

cat(sprintf("Task %d: setting %d (n=%d, p=%d, s0=%d), rep %d\n",
            task_id, setting_id, n, p, s0, rep_id))

# ================================================================
#  PARAMETERS
# ================================================================
L_outer_max <- 50
tol_outer   <- 1e-4

M_imp    <- 20                                          # number of chains
T_burn_chain <- max(20, ceiling(3 * log(n * p)))        # burn-in per chain
S_inner  <- max(5, ceiling(log(n * p)))

cat(sprintf("  Params: M=%d, T_burn=%d, S_inner=%d\n",
            M_imp, T_burn_chain, S_inner))

nfold <- 5; n_lambda <- 100; tol_lasso <- 1e-6; maxiter_lasso <- 50000

beta0 <- rep(0, p)
beta0[1] <- 1.0
if(p >= 5) beta0[2:5] <- c(1.0, 1.0, 0.5, 0.5)

n_miss_vars <- max(2, floor(p * 0.2))
miss_var_idx <- 2:min(n_miss_vars + 1, p)

# ================================================================
#  DATA GENERATION
# ================================================================
generate_data <- function(n, p, s0, rho) {
  Sigma <- rho^abs(outer(1:p, 1:p, "-"))
  X_full <- pmax(pmin(mvrnorm(n, rep(0, p), Sigma), 5.0), -5.0)
  lp <- X_full %*% beta0
  T_surv <- rexp(n, rate = exp(lp))
  C_surv <- rexp(n, rate = 0.1)
  Y <- pmin(T_surv, C_surv)
  delta <- as.numeric(T_surv <= C_surv)
  R <- matrix(1L, n, p)
  for(j in miss_var_idx) {
    logit <- -1.5 + 0.8 * X_full[, 1]
    R[, j] <- rbinom(n, 1, pmin(pmax(1/(1+exp(-logit)), 0.80), 0.95))
  }
  X_obs <- X_full; X_obs[R == 0] <- NA
  list(Y=as.numeric(Y), delta=as.numeric(delta),
       X_full=X_full, X_obs=X_obs, R=R)
}

get_breslow <- function(Y, delta, X, beta)
  breslow_estimator_cpp(as.numeric(Y), as.numeric(delta), X, beta)

init_missing <- function(X_obs, R) {
  Xc <- X_obs
  for(j in 1:ncol(X_obs)) {
    mi <- which(R[,j]==0); ob <- which(R[,j]==1)
    if(length(mi)>0 && length(ob)>0)
      Xc[mi,j] <- sample(X_obs[ob,j], length(mi), replace=TRUE)
  }
  Xc
}

fit_cox_lasso <- function(X, Y, delta, lambda=NULL) {
  surv_y <- cbind(time=Y, status=delta)
  if(is.null(lambda)) {
    cvobj <- cv.glmnet(x=X, y=surv_y, family="cox", alpha=1,
                       standardize=FALSE, nfolds=nfold, nlambda=n_lambda)
    lambda <- cvobj$lambda.min
  }
  beta <- as.vector(coef(glmnet(x=X, y=surv_y, family="cox", alpha=1,
                                lambda=lambda, standardize=FALSE,
                                thresh=tol_lasso, maxit=maxiter_lasso)))
  list(beta=beta, lambda=lambda)
}

# ================================================================
#  DEBIASED LASSO 
# ================================================================

# Helper: run QP debiasing for a given multiplier on given data
debias_at_multiplier <- function(X, Y, delta, beta_L1, neg_dll, ssq, r, multiplier) {
  n_loc <- nrow(X); p_loc <- ncol(X)
  gamma_n <- multiplier * sqrt(log(p_loc) / n_loc)
  my_pos <- which(r$values > 0); my_rank <- length(my_pos)
  beta_db <- rep(NA, p_loc); se_db <- rep(NA, p_loc)
  Theta <- matrix(NA, p_loc, p_loc)
  if(my_rank > 0) {
    Dmat <- diag(r$values[my_pos], nrow=my_rank)
    dvec <- rep(0, my_rank)
    Amat <- t(rbind(-r$vectors[,my_pos] %*% Dmat,
                    r$vectors[,my_pos] %*% Dmat))
    for(j in 1:p_loc) {
      e_j <- rep(0, p_loc); e_j[j] <- 1
      bvec <- c(-e_j, e_j) - gamma_n * rep(1, 2*p_loc)
      res <- tryCatch(solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec),
                      error=function(e) NULL)
      if(is.null(res)) { m <- solve(ssq + gamma_n*diag(p_loc), e_j)
      } else {
        m <- as.vector(r$vectors[,my_pos] %*% res$solution)
        if(my_rank < p_loc)
          m <- m + as.vector(r$vectors[,-my_pos,drop=FALSE] %*%
                               rep(0, p_loc-my_rank))
      }
      beta_db[j] <- beta_L1[j] - sum(m * neg_dll)
      se_db[j] <- sqrt(max(m[j], 1e-12) / n_loc)
      Theta[j,] <- m
    }
  } else {
    Theta <- solve(ssq + 0.1*diag(p_loc))
    beta_db <- as.vector(beta_L1 - Theta %*% neg_dll)
    se_db <- sqrt(pmax(diag(Theta), 1e-12) / n_loc)
  }
  list(beta_db=beta_db, se=se_db, Theta=Theta)
}

# CV procedure for tuning a 
tune_qp_multiplier <- function(X, Y, delta, n_cv=5, alpha_cv=0.1) {
  n_loc <- nrow(X); p_loc <- ncol(X)
  surv_y <- cbind(time=Y, status=delta)
  
  a_candidates <- exp(seq(from=log(0.01), to=log(1.5), length.out=2))
  cv_idx <- rep(1:n_cv, length.out=n_loc)
  cv_idx <- sample(cv_idx)
  
  all_cvpl <- rep(0, length(a_candidates))
  
  for(aa in seq_along(a_candidates)) {
    a_val <- a_candidates[aa]
    cvpl <- 0
    
    for(k in 1:n_cv) {
      test_idx <- which(cv_idx == k)
      train_x <- X[-test_idx, , drop=FALSE]
      test_x <- X[test_idx, , drop=FALSE]
      train_y <- Y[-test_idx]; test_y <- Y[test_idx]
      train_d <- delta[-test_idx]; test_d <- delta[test_idx]
      train_surv <- cbind(time=train_y, status=train_d)
      n_train <- nrow(train_x)
      
      # Lasso on training data
      cv_train <- tryCatch(
        cv.glmnet(x=train_x, y=train_surv, family="cox", alpha=1,
                  standardize=FALSE, nfolds=min(nfold, n_train-1), nlambda=n_lambda),
        error=function(e) NULL)
      if(is.null(cv_train)) next
      beta_train <- as.vector(coef(glmnet(x=train_x, y=train_surv, family="cox",
                                          alpha=1, lambda=cv_train$lambda.min,
                                          standardize=FALSE, thresh=tol_lasso,
                                          maxit=maxiter_lasso)))
      
      # Score and Sigma on training data
      neg_ll_t <- 0; neg_dll_t <- rep(0, p_loc)
      neg_ddll_t <- matrix(0, p_loc, p_loc); ssq_t <- matrix(0, p_loc, p_loc)
      neg_loglik_functions_cpp_ext(neg_ll_t, neg_dll_t, neg_ddll_t, ssq_t,
                                   train_x, train_y, train_d, beta_train)
      r_t <- eigen(ssq_t); r_t$values[r_t$values <= 1e-14] <- 0
      
      # Debias on training data at this a
      db_t <- debias_at_multiplier(train_x, train_y, train_d,
                                   beta_train, neg_dll_t, ssq_t, r_t, a_val)
      
      # Hard-threshold for denoising
      pval <- 2 * pnorm(abs(db_t$beta_db / db_t$se), lower.tail=FALSE)
      beta_ht <- db_t$beta_db * as.numeric(pval < (alpha_cv / p_loc))
      
      # Evaluate log partial likelihood on test data
      cvpl <- cvpl + loglik_cpp_ext(X=test_x, time=test_y,
                                    delta=test_d, beta=beta_ht)
    }
    all_cvpl[aa] <- cvpl
  }
  
  cat(sprintf("    CV scores for a: %s\n", 
              paste(sprintf("%.2f:%.4f", a_candidates, all_cvpl), collapse=", ")))
  
  # Pick a that maximizes CV log partial likelihood
  range_cvpl <- max(all_cvpl) - min(all_cvpl)
  if(range_cvpl < 1.1 * abs(mean(all_cvpl))) {
    best_a <- 0.5  # CV is flat, use default
  } else {
    best_a <- a_candidates[which.max(all_cvpl)]
  }
  best_a
}

# Main debiased lasso function: tune a, then debias
debiased_lasso_cox <- function(X, Y, delta, multiplier=NULL) {
  n_loc <- nrow(X); p_loc <- ncol(X)
  surv_y <- cbind(time=Y, status=delta)
  
  # Tune a if not provided
  if(is.null(multiplier)) {
    multiplier <- tune_qp_multiplier(X, Y, delta)
  }
  
  # Lasso on full data
  cvobj <- cv.glmnet(x=X, y=surv_y, family="cox", alpha=1,
                     standardize=FALSE, nfolds=nfold, nlambda=n_lambda)
  beta_L1 <- as.vector(coef(glmnet(x=X, y=surv_y, family="cox", alpha=1,
                                   lambda=cvobj$lambda.min, standardize=FALSE,
                                   thresh=tol_lasso, maxit=maxiter_lasso)))
  
  # Score and Sigma on full data
  neg_ll <- 0; neg_dll <- rep(0, p_loc)
  neg_ddll <- matrix(0, p_loc, p_loc); ssq <- matrix(0, p_loc, p_loc)
  neg_loglik_functions_cpp_ext(neg_ll, neg_dll, neg_ddll, ssq, X, Y, delta, beta_L1)
  r <- eigen(ssq); r$values[r$values <= 1e-14] <- 0
  
  # Debias at chosen multiplier
  db <- debias_at_multiplier(X, Y, delta, beta_L1, neg_dll, ssq, r, multiplier)
  
  list(beta_db=db$beta_db, Theta=db$Theta, beta_L1=beta_L1, se=db$se,
       lambda=cvobj$lambda.min, multiplier=multiplier)
}

# ================================================================
#  RIDGE LAMBDA TUNING
#  λ_ridge = b * p * log(p) / n
#  Tune b via K-fold CV on prediction MSE of ridge imputation model.
#  Uses first missing variable as representative.
# ================================================================
tune_ridge_multiplier <- function(Xc, R, miss_vars, n_loc, p_loc, n_cv=5) {
  j <- miss_vars[1] + 1L
  obs_idx <- which(R[, j] == 1)
  n_obs <- length(obs_idx)
  pm1 <- p_loc - 1
  
  y_obs <- Xc[obs_idx, j]
  Z_obs <- Xc[obs_idx, -j, drop=FALSE]
  
  base_rate <- p_loc * log(max(p_loc, 2)) / n_loc
  b_upper <- max(2, 10 * sqrt(p_loc / n_loc))
  b_candidates <- exp(seq(from=log(0.01), to=log(b_upper), length.out=2))
  
  
  cv_idx <- rep(1:n_cv, length.out=n_obs)
  cv_idx <- sample(cv_idx)
  
  all_cv_mse <- rep(0, length(b_candidates))
  
  for(bb in seq_along(b_candidates)) {
    b_val <- b_candidates[bb]
    lam_c <- b_val * base_rate
    cv_mse <- 0
    
    for(k in 1:n_cv) {
      test_ii <- which(cv_idx == k)
      train_ii <- which(cv_idx != k)
      
      y_train <- y_obs[train_ii]; y_test <- y_obs[test_ii]
      Z_train <- Z_obs[train_ii, , drop=FALSE]; Z_test <- Z_obs[test_ii, , drop=FALSE]
      n_train <- length(train_ii)
      
      # Center using training data
      y_mean <- mean(y_train)
      Z_mean <- colMeans(Z_train)
      yc_train <- y_train - y_mean
      Zc_train <- sweep(Z_train, 2, Z_mean)
      Zc_test <- sweep(Z_test, 2, Z_mean)
      
      # Ridge fit on training
      ZtZ <- crossprod(Zc_train) / n_train
      diag(ZtZ) <- diag(ZtZ) + lam_c
      Zty <- crossprod(Zc_train, yc_train) / n_train
      gamma_hat <- tryCatch(solve(ZtZ, Zty), error=function(e) rep(0, pm1))
      
      # Predict on test
      yhat_test <- y_mean + as.vector(Zc_test %*% gamma_hat)
      cv_mse <- cv_mse + sum((y_test - yhat_test)^2)
    }
    all_cv_mse[bb] <- cv_mse / n_obs
  }
  
  cat(sprintf("    CV scores for b: %s\n", 
              paste(sprintf("%.2f:%.4f", b_candidates, all_cv_mse), collapse=", ")))
  range_cv <- max(all_cv_mse) - min(all_cv_mse)
  if(range_cv < 1.1 * abs(mean(all_cv_mse))) {
    best_b <- 0.005  # CV is flat, use default
  } else {
    best_b <- b_candidates[which.min(all_cv_mse)]
  }
  best_b
}

# ================================================================
#  Helper: randomize sweep order
# ================================================================
permute_sweep <- function(miss_vars, cv_lam) {
  perm <- sample.int(length(miss_vars))
  list(mv = miss_vars[perm], cl = cv_lam[perm])
}

# ================================================================
#  METHOD 1: ORACLE
# ================================================================
oracle_full <- function(dat) {
  r <- debiased_lasso_cox(dat$X_full, dat$Y, dat$delta)
  sel <- as.numeric(r$beta_L1 != 0)          # lasso active set on full data
  list(beta_hat=r$beta_db, se=r$se,
       ci_lo=r$beta_db - 1.96*r$se, ci_hi=r$beta_db + 1.96*r$se,
       sel_freq=sel, selected=sel, beta_L1=r$beta_L1)
}

# ================================================================
#  METHOD 2: SMC-DBL
# ================================================================
iro_ridge_rubin <- function(dat) {
  Y <- dat$Y; delta <- dat$delta; X_obs <- dat$X_obs; R <- dat$R
  n_loc <- nrow(X_obs); p_loc <- ncol(X_obs)
  miss_vars <- as.integer(which(colSums(R == 0) > 0) - 1L)
  lambda_j <- 0.5 * sqrt(log(max(p_loc, 2)) / n_loc)
  R_int <- matrix(as.integer(R), n_loc, p_loc)
  n_mv <- length(miss_vars)
  
  # ---- Initialize ----
  Xc <- init_missing(X_obs, R)
  cx <- fit_cox_lasso(Xc, Y, delta)
  beta_cur <- cx$beta; lambda_n <- cx$lambda
  bres <- get_breslow(Y, delta, Xc, beta_cur)
  
  # Ridge lambda: b * p * log(p) / n, use simple default for Phase 1
  base_rate <- p_loc * log(max(p_loc, 2)) / n_loc
  ridge_b_init <- 0.1
  cv_lam_init <- rep(ridge_b_init * base_rate, n_mv)
  
  # ============================================================
  #  PHASE 1: IRO Convergence
  # ============================================================
  t_phase1 <- proc.time()
  for(l in 1:L_outer_max) {
    for(s in 1:S_inner) {
      ps <- permute_sweep(miss_vars, cv_lam_init)
      smc_fcs_sweep_ridge_bayes_cpp(Xc, Y, delta, R_int, beta_cur,
                                    bres$times, bres$values,
                                    lambda_j, ps$mv, ps$cl)
    }
    
    beta_new <- as.vector(coef(glmnet(
      x=Xc, y=cbind(time=Y, status=delta), family="cox", alpha=1,
      lambda=lambda_n, standardize=FALSE,
      thresh=tol_lasso, maxit=maxiter_lasso)))
    bres <- get_breslow(Y, delta, Xc, beta_new)
    
    if(sum(abs(beta_new - beta_cur)) < tol_outer) {
      cat(sprintf("    IRO Phase 1 converged at iter %d\n", l))
      beta_cur <- beta_new; break
    }
    beta_cur <- beta_new
  }
  elapsed_phase1 <- (proc.time() - t_phase1)[3]
  cat(sprintf("    Phase 1 time: %.1fs\n", elapsed_phase1))
  
  # ============================================================
  #  TRANSITION: Tune both a (QP) and b (ridge) on converged data
  # ============================================================
  t_tune <- proc.time()
  
  # Tune a for γ_n = a * sqrt(log(p)/n) via Xia's hard-thresholding CV
  qp_multiplier <- tune_qp_multiplier(Xc, Y, delta)
  cat(sprintf("    QP multiplier a=%.3f (gamma_n=%.4f)\n",
              qp_multiplier, qp_multiplier * sqrt(log(p_loc) / n_loc)))
  
  # Tune b for λ_ridge = b * p * log(p) / n via K-fold CV
  ridge_b <- tune_ridge_multiplier(Xc, R, miss_vars, n_loc, p_loc)
  lambda_ridge <- ridge_b * base_rate
  cv_lam <- rep(lambda_ridge, n_mv)
  
  cat(sprintf("    Ridge: b=%.3f, lambda_ridge=%.2e (base_rate=%.2e)\n",
              ridge_b, lambda_ridge, base_rate))
  cat(sprintf("    Phase 1 converged lasso L1=%.3f\n", sum(abs(beta_cur))))
  
  elapsed_tune <- (proc.time() - t_tune)[3]
  cat(sprintf("    Tuning time: %.1fs\n", elapsed_tune))
  
  # Save Phase 1 converged state
  Xc_phase1 <- Xc + 0
  beta_phase1 <- beta_cur
  lambda_n_phase1 <- lambda_n
  
  # ============================================================
  #  PHASE 2: M independent chains (serial, with per-chain timing)
  # ============================================================
  run_one_chain <- function(chain_id) {
    Xc_m <- Xc_phase1 + 0
    beta_m <- beta_phase1
    bres_m <- get_breslow(Y, delta, Xc_m, beta_m)
    
    # Build cache ONCE for this chain
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
  
  cat(sprintf("    Phase 2: Running %d chains serially (%d iters × %d sweeps each)...\n",
              M_imp, T_burn_chain, S_inner))
  
  t_phase2 <- proc.time()
  
  set.seed(rep_id * 1000 + setting_id)
  chain_seeds <- sample.int(1e7, M_imp)
  
  chain_results <- vector("list", M_imp)
  for(m in 1:M_imp) {
    t_m <- proc.time()
    set.seed(chain_seeds[m])
    chain_results[[m]] <- run_one_chain(m)
    elapsed_m <- (proc.time() - t_m)[3]
    if(m == 1) {
      cat(sprintf("    Chain 1 time: %.2fs (projected total: %.1fs)\n",
                  elapsed_m, elapsed_m * M_imp))
    }
  }
  
  elapsed_phase2 <- (proc.time() - t_phase2)[3]
  cat(sprintf("    Phase 2 total time: %.1fs (avg %.2fs/chain)\n",
              elapsed_phase2, elapsed_phase2 / M_imp))
  
  # ============================================================
  #  PHASE 3: Debias each chain's dataset + Rubin pool
  # ============================================================
  t_phase3 <- proc.time()
  
  valid <- !sapply(chain_results, is.null)
  M_actual <- sum(valid)
  if(M_actual < 2) {
    cat("    ERROR: fewer than 2 chains succeeded\n")
    return(NULL)
  }
  chain_results <- chain_results[valid]
  
  bdb <- matrix(NA, M_actual, p_loc)
  selM <- matrix(NA, M_actual, p_loc)          # per-chain lasso active sets
  Thl <- vector("list", M_actual)
  for(m in 1:M_actual) {
    r <- debiased_lasso_cox(chain_results[[m]]$Xc, Y, delta,
                            multiplier=qp_multiplier)
    bdb[m, ] <- r$beta_db
    selM[m, ] <- as.numeric(r$beta_L1 != 0)
    Thl[[m]] <- r$Theta
  }
  
  bb <- colMeans(bdb)
  
  VW <- rep(0, p_loc)
  for(m in 1:M_actual)
    VW <- VW + pmax(diag(Thl[[m]]), 1e-12) / n_loc
  VW <- VW / M_actual
  
  VB <- apply(bdb, 2, var)
  V_total <- VW + (1 + 1/M_actual) * VB
  se <- sqrt(pmax(V_total, 1e-12))
  
  rm_ratio <- (1 + 1/M_actual) * VB / pmax(VW, 1e-12)
  df_rubin <- pmax((M_actual - 1) * (1 + 1/rm_ratio)^2, 2)
  
  elapsed_phase3 <- (proc.time() - t_phase3)[3]
  cat(sprintf("    Phase 3 time: %.1fs (%.2fs/chain for debiasing)\n",
              elapsed_phase3, elapsed_phase3 / M_actual))
  cat(sprintf("    TOTAL IRO-Ridge-Rubin: Phase1=%.1fs + Tune=%.1fs + Phase2=%.1fs + Phase3=%.1fs = %.1fs\n",
              elapsed_phase1, elapsed_tune, elapsed_phase2, elapsed_phase3,
              elapsed_phase1 + elapsed_tune + elapsed_phase2 + elapsed_phase3))
  
  # Variable selection: pool per-chain lasso active sets by majority vote
  # (variable is "selected" if it enters the lasso model in > 50% of the
  #  M imputed datasets)
  sel_freq <- colMeans(selM)
  selected <- as.numeric(sel_freq > 0.5)

  list(beta_hat=bb, se=se,
       ci_lo=bb - qt(0.975, df_rubin)*se,
       ci_hi=bb + qt(0.975, df_rubin)*se,
       sel_freq=sel_freq, selected=selected)
}

# ================================================================
#  METHOD 3: STANDARD IRO
# ================================================================
standard_iro_liang <- function(dat) {
  Y <- dat$Y; delta <- dat$delta; X_obs <- dat$X_obs; R <- dat$R
  n_loc <- nrow(X_obs); p_loc <- ncol(X_obs)
  miss_vars <- as.integer(which(colSums(R == 0) > 0) - 1L)
  lambda_j <- 0.5 * sqrt(log(max(p_loc, 2)) / n_loc)
  R_int <- matrix(as.integer(R), n_loc, p_loc)
  
  T_total_iro <- max(20, ceiling(3 * log(n * p))) + M_imp * max(10, ceiling(2 * log(n * p)))
  T_burn_iro <- max(20, ceiling(3 * log(n * p)))
  
  cat(sprintf("    Std IRO: T_total=%d, T_burn=%d, S_inner=%d\n",
              T_total_iro, T_burn_iro, S_inner))
  
  Xc <- init_missing(X_obs, R)
  cx <- fit_cox_lasso(Xc, Y, delta)
  beta_cur <- cx$beta; lambda_n <- cx$lambda
  bres <- get_breslow(Y, delta, Xc, beta_cur)
  
  beta_chain <- matrix(NA, T_total_iro, p_loc)
  for(t in 1:T_total_iro) {
    for(s in 1:S_inner)
      smc_fcs_sweep_cpp(Xc, Y, delta, R_int, beta_cur,
                        bres$times, bres$values, lambda_j, miss_vars)
    if(t %% 10 == 1) {
      cx <- fit_cox_lasso(Xc, Y, delta); lambda_n <- cx$lambda
    }
    beta_cur <- as.vector(coef(glmnet(
      x=Xc, y=cbind(time=Y, status=delta), family="cox", alpha=1,
      lambda=lambda_n, standardize=FALSE,
      thresh=tol_lasso, maxit=maxiter_lasso)))
    bres <- get_breslow(Y, delta, Xc, beta_cur)
    if(t > T_burn_iro) beta_chain[t - T_burn_iro, ] <- beta_cur
  }
  
  n_collected <- T_total_iro - T_burn_iro
  beta_chain <- beta_chain[1:n_collected, , drop=FALSE]
  beta_bar <- colMeans(beta_chain)
  
  neg_ll <- 0; neg_dll <- rep(0, p_loc)
  neg_ddll <- matrix(0, p_loc, p_loc); ssq <- matrix(0, p_loc, p_loc)
  neg_loglik_functions_cpp_ext(neg_ll, neg_dll, neg_ddll, ssq,
                               Xc, Y, delta, beta_cur)
  hess_inv <- tryCatch(solve(neg_ddll + 0.01*diag(p_loc)),
                       error=function(e) solve(ssq + 0.1*diag(p_loc)))
  VW <- pmax(diag(hess_inv), 1e-12) / n_loc
  VB <- apply(beta_chain, 2, var)
  V_total <- VW + VB / n_collected
  se <- sqrt(pmax(V_total, 1e-12))
  
  # Variable selection: frequency with which each variable is active in the
  # post-burn-in lasso draws; majority vote (> 50%) defines the final model
  sel_freq <- colMeans(beta_chain != 0)
  selected <- as.numeric(sel_freq > 0.5)
  
  list(beta_hat=beta_bar, se=se,
       ci_lo=beta_bar - 1.96*se, ci_hi=beta_bar + 1.96*se,
       sel_freq=sel_freq, selected=selected)
}

# ================================================================
#  METHOD 4: STANDARD SMC-FCS + DEBIAS + RUBIN
# ================================================================
standard_smc_fcs <- function(dat) {
  Y <- dat$Y; delta <- dat$delta; X_obs <- dat$X_obs; R <- dat$R
  n_loc <- nrow(X_obs); p_loc <- ncol(X_obs)
  
  df <- as.data.frame(X_obs)
  colnames(df) <- paste0("X", 1:p_loc)
  df$time <- Y; df$status <- as.integer(delta)
  
  xvars <- paste0("X", 1:p_loc)
  sm_formula <- paste0("Surv(time, status) ~ ", paste(xvars, collapse=" + "))
  
  method_vec <- rep("", ncol(df))
  names(method_vec) <- colnames(df)
  miss_cols <- which(colSums(R == 0) > 0)
  for(j in miss_cols) method_vec[paste0("X", j)] <- "norm"
  
  imps <- smcfcs::smcfcs(
    originaldata = df, smtype = "coxph", smformula = sm_formula,
    method = method_vec, m = M_imp, numit = 20,
    rjlimit = 1000, noisy = FALSE
  )
  
  bdb <- matrix(NA, M_imp, p_loc)
  selM <- matrix(NA, M_imp, p_loc)             # per-imputation lasso active sets
  Thl <- vector("list", M_imp)
  for(m in 1:M_imp) {
    Xm <- as.matrix(imps$impDatasets[[m]][, xvars])
    r <- debiased_lasso_cox(Xm, Y, delta)
    bdb[m, ] <- r$beta_db; Thl[[m]] <- r$Theta
    selM[m, ] <- as.numeric(r$beta_L1 != 0)
  }
  
  bb <- colMeans(bdb)
  VW <- rep(0, p_loc)
  for(m in 1:M_imp) VW <- VW + pmax(diag(Thl[[m]]), 1e-12) / n_loc
  VW <- VW / M_imp
  VB <- apply(bdb, 2, var)
  V_total <- VW + (1 + 1/M_imp) * VB
  se <- sqrt(pmax(V_total, 1e-12))
  rm_ratio <- (1 + 1/M_imp) * VB / pmax(VW, 1e-12)
  df_rubin <- pmax((M_imp - 1) * (1 + 1/rm_ratio)^2, 2)
  
  sel_freq <- colMeans(selM)
  selected <- as.numeric(sel_freq > 0.5)
  
  list(beta_hat=bb, se=se,
       ci_lo=bb - qt(0.975, df_rubin)*se,
       ci_hi=bb + qt(0.975, df_rubin)*se,
       sel_freq=sel_freq, selected=selected)
}

# ================================================================
#  RUN
# ================================================================
set.seed(rep_id * 137 + 7)
dat <- generate_data(n, p, s0, rho)

methods <- c("oracle", "iro_ridge_rubin", "standard_iro")
method_fns <- list(
  oracle          = oracle_full,
  iro_ridge_rubin = iro_ridge_rubin,
  standard_iro    = standard_iro_liang
)

res_all <- list()
for(method in methods) {
  cat(sprintf("  Running %s...\n", method))
  t0 <- proc.time()
  res <- tryCatch(method_fns[[method]](dat), error=function(e) {
    cat(sprintf("    FAILED: %s\n", e$message)); NULL })
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("    Done (%.1fs)\n", elapsed))
  if(!is.null(res)) {
    res_all[[method]] <- list(beta_hat=res$beta_hat, se=res$se,
                              ci_lo=res$ci_lo, ci_hi=res$ci_hi,
                              sel_freq=res$sel_freq, selected=res$selected)
    # Per-replicate variable-screening summary (penalized selection):
    #   MS = selected model size, TP = true positives, FP = false positives
    if(!is.null(res$selected)) {
      S_true <- which(beta0 != 0)
      sel_idx <- which(res$selected == 1)
      cat(sprintf("    Selection: MS=%d, TP=%d, FP=%d\n",
                  length(sel_idx),
                  sum(sel_idx %in% S_true),
                  sum(!(sel_idx %in% S_true))))
    }
  }
}

dir.create("results", showWarnings = FALSE)
outfile <- sprintf("results/res_s%02d_r%03d.RData", setting_id, rep_id)
save(res_all, n, p, s0, rho, beta0, setting_id, rep_id, file = outfile)
cat(sprintf("Saved to %s\n", outfile))