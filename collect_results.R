###############################################################################
#  collect_results.R
###############################################################################
settings <- expand.grid(n = c(500,1000,2000), p = c(20,50,100,200))
settings$s0 <- 5; n_settings <- nrow(settings); n_reps <- 100
methods <- c("oracle","iro_ridge_rubin","standard_iro","std_smc_fcs","mean_imp")
labels  <- c("Oracle","SMC-DBL","Standard IRO (Liang)","Std SMC-FCS","Mean Imp + DL")

table1 <- data.frame()
table2 <- data.frame()
table3 <- data.frame()   # variable-screening: model size, TP, FP
percoef <- data.frame()

for(sid in 1:n_settings) {
  nn <- settings$n[sid]; pp <- settings$p[sid]
  beta0 <- rep(0, pp); beta0[1] <- 1
  if(pp >= 5) beta0[2:5] <- c(1, 1, 0.5, 0.5)
  S <- which(beta0 != 0); Z <- which(beta0 == 0)

  loaded <- 0; store <- list()
  for(r in 1:n_reps) {
    f <- sprintf("results/res_s%02d_r%03d.RData", sid, r)
    if(!file.exists(f)) next
    load(f); loaded <- loaded + 1; store[[loaded]] <- res_all
  }
  if(loaded < 5) next

  for(mi in seq_along(methods)) {
    mn <- methods[mi]; ml <- labels[mi]
    bh <- se_m <- cv <- matrix(NA, loaded, pp)
    typeI_z <- matrix(NA, loaded, pp)
    selmat <- matrix(NA, loaded, pp)   # 0/1 penalized selection per replicate
    nv <- 0
    for(r in 1:loaded) {
      if(is.null(store[[r]][[mn]])) next
      nv <- nv + 1; res <- store[[r]][[mn]]
      bh[nv, ] <- res$beta_hat; se_m[nv, ] <- res$se
      cv[nv, ] <- as.numeric(res$ci_lo <= beta0 & beta0 <= res$ci_hi)
      # Type I error indicator: reject H0: beta_k = 0 at level 0.05
      typeI_z[nv, ] <- as.numeric(!(res$ci_lo <= 0 & 0 <= res$ci_hi))
      # Penalized selection indicator (NA for legacy result files that
      # predate the variable-screening addition)
      if(!is.null(res$selected)) selmat[nv, ] <- res$selected
    }
    if(nv < 5) next
    bh <- bh[1:nv, , drop = FALSE]
    se_m <- se_m[1:nv, , drop = FALSE]
    cv <- cv[1:nv, , drop = FALSE]
    typeI_z <- typeI_z[1:nv, , drop = FALSE]
    selmat <- selmat[1:nv, , drop = FALSE]

    bias_k <- colMeans(bh) - beta0
    rmse_k <- sqrt(colMeans((bh - rep(1, nv) %o% beta0)^2))
    emp_k <- apply(bh, 2, sd)
    avg_k <- colMeans(se_m)
    cov_k <- colMeans(cv)
    typeI_k <- colMeans(typeI_z)

    # Active-set averages
    abs_bias_S <- mean(abs(bias_k[S]))
    rmse_S     <- mean(rmse_k[S])
    emp_S      <- mean(emp_k[S])
    avg_S      <- mean(avg_k[S])
    ratio_S    <- avg_S / max(emp_S, 1e-12)
    cov_S      <- mean(cov_k[S])

    # Null-set averages (only if there are any nulls)
    if(length(Z) > 0) {
      emp_Z   <- mean(emp_k[Z])
      avg_Z   <- mean(avg_k[Z])
      ratio_Z <- avg_Z / max(emp_Z, 1e-12)
      cov_Z   <- mean(cov_k[Z])
      typeI   <- mean(typeI_k[Z])
    } else {
      emp_Z <- avg_Z <- ratio_Z <- cov_Z <- typeI <- NA
    }

    table1 <- rbind(table1, data.frame(
      n = nn, p = pp, method = ml, n_valid = nv,
      abs_bias_S = abs_bias_S, rmse_S = rmse_S,
      cov_S = cov_S, cov_Z = cov_Z, typeI_Z = typeI
    ))

    table2 <- rbind(table2, data.frame(
      n = nn, p = pp, method = ml, n_valid = nv,
      emp_S = emp_S, avg_S = avg_S, ratio_S = ratio_S,
      emp_Z = emp_Z, avg_Z = avg_Z, ratio_Z = ratio_Z
    ))

    # ---- Table 3: variable-screening performance of the penalized fit ----
    sel_ok <- which(rowSums(is.na(selmat)) == 0)   # replicates with selection info
    if(length(sel_ok) >= 5) {
      sm <- selmat[sel_ok, , drop = FALSE]
      ms_r <- rowSums(sm)                                   # model size per rep
      tp_r <- rowSums(sm[, S, drop = FALSE])                # true positives
      fp_r <- if(length(Z) > 0) rowSums(sm[, Z, drop = FALSE]) else rep(0, nrow(sm))
      sure_r  <- as.numeric(tp_r == length(S))              # S ⊆ Ŝ
      exact_r <- as.numeric(tp_r == length(S) & fp_r == 0)  # Ŝ = S

      table3 <- rbind(table3, data.frame(
        n = nn, p = pp, method = ml, n_valid = length(sel_ok),
        MS = mean(ms_r), MS_sd = sd(ms_r),
        TP = mean(tp_r), TP_sd = sd(tp_r),
        FP = mean(fp_r), FP_sd = sd(fp_r),
        sure_screen = mean(sure_r),   # P(all truly non-zero retained)
        exact_sel   = mean(exact_r)   # P(selected set equals true set)
      ))
    }

    for(k in 1:pp) {
      percoef <- rbind(percoef, data.frame(
        n = nn, p = pp, method = ml, k = k, beta0 = beta0[k],
        bias = bias_k[k], rmse = rmse_k[k], emp = emp_k[k],
        avg = avg_k[k], cov = cov_k[k]
      ))
    }
  }
}

dir.create("simout", showWarnings = FALSE)
write.csv(table1,  "simout/table1_headline.csv",      row.names = FALSE)
write.csv(table2,  "simout/table2_variance_calib.csv", row.names = FALSE)
write.csv(table3,  "simout/table3_selection.csv",      row.names = FALSE)
write.csv(percoef, "simout/percoef_long.csv",          row.names = FALSE)

cat("Wrote simout/table1_headline.csv, simout/table2_variance_calib.csv,",
    "simout/table3_selection.csv, simout/percoef_long.csv\n")
