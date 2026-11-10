###############################################################################
#  bench_complete.R
#
#  Standalone runner for a baseline comparison: column-wise mean imputation
#  followed by a single debiased-lasso fit on the imputed design matrix. 
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])

library(Rcpp); library(RcppArmadillo); library(glmnet)
library(quadprog); library(MASS); library(survival)
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
#  PARAMETERS  (only those needed for the lasso/debias call)
# ================================================================
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

# ================================================================
#  DEBIASED LASSO 
# ================================================================
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
      mSm <- as.numeric(t(m) %*% ssq %*% m)
      mSm <- max(mSm, 1e-12)
      se_db[j] <- sqrt(mSm / n_loc)
      Theta[j,] <- m
    }
  } else {
    Theta <- solve(ssq + 0.1*diag(p_loc))
    beta_db <- as.vector(beta_L1 - Theta %*% neg_dll)
    var_diag <- diag(Theta %*% ssq %*% t(Theta)) / n_loc
    se_db <- sqrt(pmax(var_diag, 1e-12))
  }
  list(beta_db=beta_db, se=se_db, Theta=Theta)
}

debiased_lasso_cox <- function(X, Y, delta, multiplier=0.5) {
  n_loc <- nrow(X); p_loc <- ncol(X)
  surv_y <- cbind(time=Y, status=delta)

  cvobj <- cv.glmnet(x=X, y=surv_y, family="cox", alpha=1,
                     standardize=FALSE, nfolds=nfold, nlambda=n_lambda)
  beta_L1 <- as.vector(coef(glmnet(x=X, y=surv_y, family="cox", alpha=1,
                                    lambda=cvobj$lambda.min, standardize=FALSE,
                                    thresh=tol_lasso, maxit=maxiter_lasso)))

  neg_ll <- 0; neg_dll <- rep(0, p_loc)
  neg_ddll <- matrix(0, p_loc, p_loc); ssq <- matrix(0, p_loc, p_loc)
  neg_loglik_functions_cpp_ext(neg_ll, neg_dll, neg_ddll, ssq, X, Y, delta, beta_L1)
  r <- eigen(ssq); r$values[r$values <= 1e-14] <- 0

  db <- debias_at_multiplier(X, Y, delta, beta_L1, neg_dll, ssq, r, multiplier)

  list(beta_db=db$beta_db, Theta=db$Theta, beta_L1=beta_L1, se=db$se,
       ssq=ssq, lambda=cvobj$lambda.min, multiplier=multiplier)
}

# ================================================================
#  EMPTY RESULT
#
#  Always returns vectors of length p.
#  This guarantees identical output dimensions even if there are
#  too few complete cases to fit the model.
# ================================================================
empty_result <- function(p) {
  
  list(
    beta_hat = rep(NA_real_, p),
    se       = rep(NA_real_, p),
    ci_lo    = rep(NA_real_, p),
    ci_hi    = rep(NA_real_, p),
    sel_freq = rep(NA_real_, p),
    selected = rep(NA_real_, p)
  )
}

# ================================================================
#  METHOD: COMPLETE CASE + DEBIASED LASSO
# ================================================================
complete_case_debias <- function(dat) {
  
  # Subject is retained only if ALL p covariates are observed
  cc_idx <- complete.cases(dat$X_obs)
  
  n_cc <- sum(cc_idx)
  
  cat(sprintf(
    "    Complete cases: %d / %d (%.2f%%)\n",
    n_cc,
    nrow(dat$X_obs),
    100 * n_cc / nrow(dat$X_obs)
  ))
  
  if(n_cc < nfold) {
    
    cat(sprintf(
      "    Too few complete cases (%d < %d); returning NA vectors.\n",
      n_cc, nfold
    ))
    
    return(empty_result(ncol(dat$X_obs)))
  }
  
  X_cc <- dat$X_obs[
    cc_idx,
    ,
    drop = FALSE
  ]
  
  Y_cc <- dat$Y[cc_idx]
  delta_cc <- dat$delta[cc_idx]
  
  # Need observed events for Cox regression
  if(sum(delta_cc) < 2) {
    
    cat(sprintf(
      "    Too few events among complete cases (%d); returning NA vectors.\n",
      sum(delta_cc)
    ))
    
    return(empty_result(ncol(dat$X_obs)))
  }
  
  r <- tryCatch(
    debiased_lasso_cox(
      X_cc,
      Y_cc,
      delta_cc
    ),
    error = function(e) {
      
      cat(sprintf(
        "    Complete-case debiased Lasso failed: %s\n",
        e$message
      ))
      
      NULL
    }
  )
  
  if(is.null(r)) {
    return(empty_result(ncol(dat$X_obs)))
  }
  sel <- as.numeric(r$beta_L1 != 0)
  list(
    beta_hat = r$beta_db,
    se       = r$se,
    ci_lo    = r$beta_db - 1.96 * r$se,
    ci_hi    = r$beta_db + 1.96 * r$se,
    sel_freq=sel, selected=sel
  )
}


# ================================================================
#  RUN
# ================================================================
set.seed(rep_id * 137 + 7)
dat <- generate_data(n, p, s0, rho)

cat("  Running complete...\n")
t0 <- proc.time()
res <- tryCatch(complete_case_debias(dat),
                error = function(e) {
                  cat(sprintf("    FAILED: %s\n", e$message)); NULL
                })
elapsed <- (proc.time() - t0)[3]
cat(sprintf("    Done (%.1fs)\n", elapsed))

# ================================================================
#  SAVE
# ================================================================
dir.create("results", showWarnings = FALSE)
outfile <- sprintf("results/res_s%02d_r%03d.RData", setting_id, rep_id)

if(file.exists(outfile)) {
  e <- new.env()
  load(outfile, envir = e)
  res_all <- if(exists("res_all", envir = e)) e$res_all else list()
} else {
  res_all <- list()
}

if(!is.null(res)) {
  res_all[["complete"]] <- list(beta_hat = res$beta_hat,
                                se       = res$se,
                                ci_lo    = res$ci_lo,
                                ci_hi    = res$ci_hi,
                                sel_freq=res$sel_freq,
                                selected=res$selected)
}

save(res_all, n, p, s0, rho, beta0, setting_id, rep_id, file = outfile)
cat(sprintf("Saved (complete appended) to %s\n", outfile))
