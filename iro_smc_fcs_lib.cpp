// iro_smc_fcs_lib.cpp
// C++ backend for IRO-SMC-FCS with Bayesian ridge imputation
//
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// ================================================================
//  Utilities
// ================================================================
static double sign_f(double a) {
  if(a == 0.0) return 0.0;
  return (a > 0.0 ? 1.0 : -1.0);
}
static double abs_f(double a) {
  return (a > 0.0 ? a : -a);
}

// ================================================================
//  Sort-based cumulative sum helpers for Cox likelihood
//  Sort observations by DECREASING time.  Walking forward through
//  this order, the risk set at event time t_i is the running sum
//  of exp(Xβ) up to and including all subjects with time >= t_i.
// ================================================================

// ================================================================
//  neg_loglik_functions (internal) — O(n·p² + n·log(n))
// ================================================================
static void neg_loglik_functions_cpp_internal(
    double& neg_loglik, arma::colvec& neg_dloglik, arma::mat& neg_ddloglik,
    arma::mat& X, arma::colvec& time, arma::colvec& delta, arma::colvec& beta) {

  int p = X.n_cols, n = X.n_rows;
  arma::colvec xbeta = X * beta;
  arma::colvec exp_xbeta = exp(xbeta);

  // Sort by DECREASING time
  arma::uvec ord = sort_index(time, "descend");

  // Running cumulative sums over the risk set
  double S0 = 0.0;
  arma::colvec S1(p, fill::zeros);
  arma::mat S2(p, p, fill::zeros);

  neg_loglik = 0.0; neg_dloglik.fill(0); neg_ddloglik.fill(0);
  double loglik_acc = arma::as_scalar(xbeta.t() * delta);

  // Walk in decreasing time order; accumulate risk set
  // When we encounter an event, the running sums are exactly the risk set sums
  // Handle ties: accumulate all subjects at same time, then process events
  int idx = 0;
  while(idx < n) {
    double t_cur = time(ord(idx));
    // Accumulate all subjects at this time point into risk set
    int idx_start = idx;
    while(idx < n && time(ord(idx)) == t_cur) {
      int i = ord(idx);
      double ew = exp_xbeta(i);
      S0 += ew;
      S1 += X.row(i).t() * ew;
      S2 += (X.row(i).t() * X.row(i)) * ew;
      idx++;
    }
    // Now process events at this time point
    for(int kk = idx_start; kk < idx; kk++) {
      int i = ord(kk);
      if(delta(i) == 0) continue;
      double s0 = S0 / (double)n;
      arma::colvec s1 = S1 / (double)n;
      arma::mat s2 = S2 / (double)n;
      loglik_acc -= log(s0);
      arma::colvec s1_over_s0 = s1 / s0;
      neg_dloglik += X.row(i).t() - s1_over_s0;
      neg_ddloglik += s2 / s0 - s1_over_s0 * s1_over_s0.t();
    }
  }
  neg_loglik = -1.0 / (double)n * loglik_acc;
  neg_dloglik = -1.0 / (double)n * neg_dloglik;
  neg_ddloglik = neg_ddloglik / (double)n;
}

// ================================================================
//  neg_loglik_functions_cpp_ext (exported, with score_sq) — FAST
// ================================================================
// [[Rcpp::export]]
void neg_loglik_functions_cpp_ext(
    double& neg_loglik, arma::colvec& neg_dloglik, arma::mat& neg_ddloglik,
    arma::mat& score_sq,
    arma::mat& X, arma::colvec& time, arma::colvec& delta, arma::colvec& beta) {

  int p = X.n_cols, n = X.n_rows;
  arma::colvec xbeta = X * beta;
  arma::colvec exp_xbeta = exp(xbeta);

  arma::uvec ord = sort_index(time, "descend");

  double S0 = 0.0;
  arma::colvec S1(p, fill::zeros);
  arma::mat S2(p, p, fill::zeros);

  neg_loglik = 0.0; neg_dloglik.fill(0); neg_ddloglik.fill(0); score_sq.fill(0);
  double loglik_acc = arma::as_scalar(xbeta.t() * delta);

  int idx = 0;
  while(idx < n) {
    double t_cur = time(ord(idx));
    int idx_start = idx;
    while(idx < n && time(ord(idx)) == t_cur) {
      int i = ord(idx);
      double ew = exp_xbeta(i);
      S0 += ew;
      S1 += X.row(i).t() * ew;
      S2 += (X.row(i).t() * X.row(i)) * ew;
      idx++;
    }
    for(int kk = idx_start; kk < idx; kk++) {
      int i = ord(kk);
      if(delta(i) == 0) continue;
      double s0 = S0 / (double)n;
      arma::colvec s1 = S1 / (double)n;
      arma::mat s2 = S2 / (double)n;
      loglik_acc -= log(s0);
      arma::colvec s1_over_s0 = s1 / s0;
      arma::colvec resid_i = X.row(i).t() - s1_over_s0;
      neg_dloglik += resid_i;
      neg_ddloglik += s2 / s0 - s1_over_s0 * s1_over_s0.t();
      score_sq += resid_i * resid_i.t();
    }
  }
  neg_loglik = -1.0 / (double)n * loglik_acc;
  neg_dloglik = -1.0 / (double)n * neg_dloglik;
  neg_ddloglik = neg_ddloglik / (double)n;
  score_sq = score_sq / (double)n;
}

// ================================================================
//  lasso_univCox_cpp_ext
// ================================================================
// [[Rcpp::export]]
arma::colvec lasso_univCox_cpp_ext(double lambda, arma::mat& X,
                                   arma::colvec& time, arma::colvec& delta,
                                   arma::colvec& beta_init,
                                   double tol = 1e-6, int maxiter = 1000) {
  int p = X.n_cols;
  arma::colvec beta_old = beta_init, beta_new = beta_init;
  double neg_loglik; arma::colvec neg_dloglik(p); arma::mat neg_ddloglik(p, p);
  for(int iter = 0; iter < maxiter; iter++) {
    neg_loglik = 0; neg_dloglik.fill(0); neg_ddloglik.fill(0);
    neg_loglik_functions_cpp_internal(neg_loglik, neg_dloglik, neg_ddloglik,
                                     X, time, delta, beta_old);
    for(int d = 0; d < p; d++) {
      double a_d = arma::as_scalar(trans(beta_new) * neg_ddloglik.col(d))
                   - beta_new(d) * neg_ddloglik(d, d)
                   - arma::as_scalar(trans(beta_old) * neg_ddloglik.col(d))
                   + neg_dloglik(d);
      beta_new(d) = -sign_f(a_d) *
        (abs_f(a_d) > lambda ? abs_f(a_d) - lambda : 0.0) /
        neg_ddloglik(d, d);
    }
    double diff = sqrt(sum(pow(beta_new - beta_old, 2)));
    beta_old = beta_new;
    if(diff < tol) break;
  }
  return beta_new;
}

// ================================================================
//  loglik_cpp_ext — FAST O(n log n + n·p)
// ================================================================
// [[Rcpp::export]]
double loglik_cpp_ext(arma::mat& X, arma::colvec& time,
                      arma::colvec& delta, arma::colvec& beta) {
  int n = X.n_rows;
  arma::colvec xbeta = X * beta;
  arma::colvec exp_xbeta = exp(xbeta);
  double loglik = arma::as_scalar(xbeta.t() * delta);

  arma::uvec ord = sort_index(time, "descend");
  double S0 = 0.0;
  int idx = 0;
  while(idx < n) {
    double t_cur = time(ord(idx));
    int idx_start = idx;
    while(idx < n && time(ord(idx)) == t_cur) {
      S0 += exp_xbeta(ord(idx));
      idx++;
    }
    for(int kk = idx_start; kk < idx; kk++) {
      if(delta(ord(kk)) == 0) continue;
      loglik -= log(S0 / (double)n);
    }
  }
  return 1.0 / (double)n * loglik;
}

// ================================================================
//  C_stat_harrell
// ================================================================
// [[Rcpp::export]]
double C_stat_harrell(arma::mat& Z, arma::colvec& beta,
                      arma::colvec& delta, arma::colvec& time) {
  int n = Z.n_rows; arma::colvec zbeta = Z * beta;
  double denom = 0.0, num = 0.0;
  for(int i = 0; i < n; i++) {
    if(delta(i) == 1.0) {
      for(int j = 0; j < n; j++) {
        if(j != i) {
          num += ((time(i) < time(j)) && (zbeta(i) > zbeta(j))) ? 1.0 : 0.0;
          denom += (time(i) < time(j)) ? 1.0 : 0.0;
        }
      }
    }
  }
  return num / denom;
}

// ================================================================
//  Breslow estimator — FAST O(n log n)
// ================================================================
// [[Rcpp::export]]
List breslow_estimator_cpp(arma::colvec& time, arma::colvec& delta,
                           arma::mat& X, arma::colvec& beta) {
  int n = X.n_rows;
  arma::colvec exp_lp = exp(clamp(X * beta, -10.0, 10.0));

  // Sort by increasing time
  arma::uvec ord = sort_index(time, "ascend");

  // Reverse cumsum: risk set sum at time t_k = sum of exp_lp for all i with time(i) >= t_k
  // Compute by walking from largest to smallest and accumulating
  arma::vec risk_at(n);
  double running_risk = 0.0;
  for(int ii = n - 1; ii >= 0; ii--) {
    running_risk += exp_lp(ord(ii));
    risk_at(ii) = running_risk;
  }

  // Collect unique event times and their Breslow increments
  arma::uvec ev_idx = find(delta == 1);
  arma::vec ev_times = sort(unique(time.elem(ev_idx)));
  int n_events = ev_times.n_elem;
  arma::vec cum_haz(n_events, fill::zeros);

  double ch = 0.0;
  int sorted_pos = 0; // position in sorted order
  for(int k = 0; k < n_events; k++) {
    double t_k = ev_times(k);
    // Advance sorted_pos to first subject with time >= t_k
    while(sorted_pos < n && time(ord(sorted_pos)) < t_k) sorted_pos++;
    double risk = (sorted_pos < n) ? risk_at(sorted_pos) : 1e-10;
    // Count events at t_k
    double d_k = 0.0;
    for(int ii = sorted_pos; ii < n && time(ord(ii)) == t_k; ii++) {
      if(delta(ord(ii)) == 1) d_k += 1.0;
    }
    ch += d_k / std::max(risk, 1e-10);
    cum_haz(k) = ch;
  }
  return List::create(_["times"] = ev_times, _["values"] = cum_haz);
}

// ================================================================
//  Individual Cox likelihood (for MH acceptance)
// ================================================================
static double cox_lik_individual(double Yi, double delta_i,
                                 const arma::rowvec& Xi,
                                 const arma::colvec& beta,
                                 const arma::vec& ch_times,
                                 const arma::vec& ch_vals) {
  double lp = arma::as_scalar(Xi * beta);
  lp = std::max(-10.0, std::min(10.0, lp));
  double Lambda0 = 0.0;
  for(int k = 0; k < (int)ch_times.n_elem; k++) {
    if(ch_times(k) <= Yi) Lambda0 = ch_vals(k);
    else break;
  }
  double log_f = delta_i * lp - exp(lp) * Lambda0;
  log_f = std::max(-20.0, std::min(20.0, log_f));
  return exp(log_f);
}

// ================================================================
//  Precompute ridge cache for Phase 2
//  Returns an R list of lists, one per missing variable.
// ================================================================
// [[Rcpp::export]]
List precompute_ridge_cache_cpp(arma::mat& Xc,
                                arma::imat& R,
                                arma::ivec& miss_vars,
                                arma::vec& cv_lam) {
  int n = Xc.n_rows, p = Xc.n_cols;
  int n_mv = (int)miss_vars.n_elem;
  List cache(n_mv);

  for(int jj = 0; jj < n_mv; jj++) {
    int j = miss_vars(jj);
    int pm1 = p - 1;

    // Partition
    std::vector<int> obs_v, miss_v;
    for(int i = 0; i < n; i++) {
      if(R(i, j) == 0) miss_v.push_back(i);
      else obs_v.push_back(i);
    }
    int n_obs = (int)obs_v.size();
    int n_miss = (int)miss_v.size();

    arma::ivec obs_idx(n_obs), miss_idx(n_miss);
    for(int i = 0; i < n_obs; i++) obs_idx(i) = obs_v[i];
    for(int i = 0; i < n_miss; i++) miss_idx(i) = miss_v[i];

    // Build Z_obs and y_obs
    arma::mat Z_obs(n_obs, pm1);
    arma::colvec y_obs(n_obs);
    for(int ii = 0; ii < n_obs; ii++) {
      int i = obs_v[ii];
      y_obs(ii) = Xc(i, j);
      int col = 0;
      for(int k = 0; k < p; k++) {
        if(k == j) continue;
        Z_obs(ii, col++) = Xc(i, k);
      }
    }

    // Center
    double y_mean = arma::mean(y_obs);
    arma::colvec yc = y_obs - y_mean;
    arma::rowvec Z_mean = arma::mean(Z_obs, 0);
    arma::mat Zc = Z_obs.each_row() - Z_mean;

    // ZtZ / n_obs + lambda * I
    double lambda_ridge = cv_lam(jj);
    arma::mat ZtZ = Zc.t() * Zc / (double)n_obs;
    ZtZ.diag() += lambda_ridge;

    // Cholesky
    arma::mat L_ZtZ;
    bool chol_ok = arma::chol(L_ZtZ, ZtZ, "lower");
    if(!chol_ok) {
      ZtZ.diag() += 0.01;
      chol_ok = arma::chol(L_ZtZ, ZtZ, "lower");
      if(!chol_ok) {
        // Last resort: use identity
        L_ZtZ = arma::eye(pm1, pm1);
      }
    }

    // Zty / n_obs (for the initial cached version)
    arma::colvec Zty = Zc.t() * yc / (double)n_obs;

    // Solve for initial gamma_hat using the Cholesky
    arma::colvec gamma_hat = arma::solve(trimatl(L_ZtZ), Zty);
    gamma_hat = arma::solve(trimatu(L_ZtZ.t()), gamma_hat);

    // Residual variance
    arma::colvec resid = y_obs - y_mean - Z_obs * gamma_hat;
    int df = std::max(n_obs - pm1 - 1, 1);
    double sigma2_hat = std::max(arma::dot(resid, resid) / (double)df, 0.001);

    cache[jj] = List::create(
      _["j"] = j,
      _["obs_idx"] = obs_idx,
      _["miss_idx"] = miss_idx,
      _["n_obs"] = n_obs,
      _["n_miss"] = n_miss,
      _["y_obs"] = y_obs,
      _["y_mean"] = y_mean,
      _["Z_obs"] = Z_obs,
      _["Z_mean"] = Z_mean,
      _["Zc"] = Zc,
      _["Zty"] = Zty,
      _["L_ZtZ"] = L_ZtZ,
      _["gamma_hat"] = gamma_hat,
      _["sigma2_hat"] = sigma2_hat,
      _["df"] = df,
      _["lambda_ridge"] = lambda_ridge
    );
  }
  return cache;
}

// ================================================================
//  Cached Bayesian ridge sweep (Phase 2)
// ================================================================
// [[Rcpp::export]]
void smc_fcs_sweep_ridge_bayes_cached_cpp(
    arma::mat& Xc,
    arma::colvec& time, arma::colvec& delta,
    arma::imat& R,
    arma::colvec& beta,
    arma::vec& ch_times, arma::vec& ch_vals,
    double lambda_j,
    arma::ivec& miss_vars,
    List& cache) {

  int n = Xc.n_rows, p = Xc.n_cols;

  for(int jj = 0; jj < (int)miss_vars.n_elem; jj++) {
    // Use miss_vars(jj) to find which cache entry
    // The cache was built in the original miss_vars order,
    // but this call may have a PERMUTED order.
    // We need to find the cache entry for variable miss_vars(jj).
    int j = miss_vars(jj);
    int pm1 = p - 1;

    // Find the cache entry for variable j
    int cache_idx = -1;
    for(int cc = 0; cc < cache.size(); cc++) {
      List c = cache[cc];
      int cj = c["j"];
      if(cj == j) { cache_idx = cc; break; }
    }
    if(cache_idx < 0) continue;

    List cj_list = cache[cache_idx];
    arma::ivec obs_idx = cj_list["obs_idx"];
    arma::ivec miss_idx = cj_list["miss_idx"];
    int n_obs = cj_list["n_obs"];
    int n_miss = cj_list["n_miss"];
    arma::colvec y_obs = cj_list["y_obs"];   // FIXED observed values
    double y_mean = cj_list["y_mean"];        // FIXED
    arma::rowvec Z_mean_cached = cj_list["Z_mean"];
    arma::mat L_ZtZ = cj_list["L_ZtZ"];      // FIXED Cholesky
    int df = cj_list["df"];

    if(n_miss == 0) continue;

    // ---- Refresh Z_obs from current Xc and recompute Zty ----
    // y_obs is fixed (observed values of col j), but Z_obs columns
    // for other variables may have changed due to imputation.
    // We recompute Zty = Zc' * yc / n_obs using current Xc.
    // This is O(n_obs * p), much cheaper than O(n_obs * p^2) for ZtZ.
    arma::colvec Zty(pm1, fill::zeros);
    arma::rowvec Z_mean(pm1, fill::zeros);

    // Compute Z_mean from current Xc (observed rows for variable j)
    for(int ii = 0; ii < n_obs; ii++) {
      int i = obs_idx(ii);
      int col = 0;
      for(int k = 0; k < p; k++) {
        if(k == j) continue;
        Z_mean(col) += Xc(i, k);
        col++;
      }
    }
    Z_mean /= (double)n_obs;

    // Compute Zty = sum_obs (Z_i - Z_mean)' * (y_i - y_mean) / n_obs
    for(int ii = 0; ii < n_obs; ii++) {
      int i = obs_idx(ii);
      double yi_c = y_obs(ii) - y_mean;
      int col = 0;
      for(int k = 0; k < p; k++) {
        if(k == j) continue;
        Zty(col) += (Xc(i, k) - Z_mean(col)) * yi_c;
        col++;
      }
    }
    Zty /= (double)n_obs;

    // Solve for gamma_hat using cached Cholesky: L L' gamma = Zty
    arma::colvec gamma_hat = arma::solve(trimatl(L_ZtZ), Zty);
    gamma_hat = arma::solve(trimatu(L_ZtZ.t()), gamma_hat);

    double alpha_hat = y_mean - arma::dot(Z_mean, gamma_hat);

    // Residual variance from observed subjects
    double ss = 0.0;
    for(int ii = 0; ii < n_obs; ii++) {
      int i = obs_idx(ii);
      double pred = alpha_hat;
      int col = 0;
      for(int k = 0; k < p; k++) {
        if(k == j) continue;
        pred += Xc(i, k) * gamma_hat(col);
        col++;
      }
      double r_i = y_obs(ii) - pred;
      ss += r_i * r_i;
    }
    double sigma2_hat = std::max(ss / (double)df, 0.001);

    // ---- BAYESIAN POSTERIOR DRAW ----
    // Draw sigma2 ~ df * sigma2_hat / chi2(df)
    double chi2_draw = R::rchisq((double)df);
    double sigma2_draw = (double)df * sigma2_hat / std::max(chi2_draw, 1e-10);
    sigma2_draw = std::max(sigma2_draw, 0.001);

    // Draw gamma ~ N(gamma_hat, sigma2_draw * (ZtZ)^{-1} / n_obs)
    // = gamma_hat + sqrt(sigma2_draw/n_obs) * L_ZtZ^{-T} * z
    // where L_ZtZ * L_ZtZ' = ZtZ, so L_ZtZ^{-T} * z has cov ZtZ^{-1}
    double scale = sqrt(sigma2_draw / (double)n_obs);
    arma::colvec z_vec(pm1);
    for(int dd = 0; dd < pm1; dd++) z_vec(dd) = R::rnorm(0.0, 1.0);
    arma::colvec gamma_draw = gamma_hat +
      scale * arma::solve(trimatu(L_ZtZ.t()), z_vec);

    double alpha_draw = y_mean - arma::dot(Z_mean, gamma_draw);
    double sigma_draw = sqrt(sigma2_draw);

    // ---- MH STEPS for each missing subject ----
    for(int ii = 0; ii < n_miss; ii++) {
      int i = miss_idx(ii);

      // Compute mu = alpha_draw + Z_i' * gamma_draw
      double mu = alpha_draw;
      int col = 0;
      for(int k = 0; k < p; k++) {
        if(k == j) continue;
        mu += Xc(i, k) * gamma_draw(col);
        col++;
      }

      double x_star = R::rnorm(mu, sigma_draw);
      x_star = std::max(-5.0, std::min(5.0, x_star));
      double x_curr = Xc(i, j);

      arma::rowvec Xi_star = Xc.row(i); Xi_star(j) = x_star;
      arma::rowvec Xi_curr = Xc.row(i); Xi_curr(j) = x_curr;

      double f_star = cox_lik_individual(time(i), delta(i), Xi_star, beta,
                                         ch_times, ch_vals);
      double f_curr = cox_lik_individual(time(i), delta(i), Xi_curr, beta,
                                         ch_times, ch_vals);

      double ratio = std::min(1.0, f_star / std::max(f_curr, 1e-30));
      if(R::runif(0, 1) < ratio) {
        Xc(i, j) = x_star;
      }
    }
  }
}

// ================================================================
//  Original Bayesian ridge sweep (kept for Phase 1 / backward compat)
// ================================================================
// [[Rcpp::export]]
void smc_fcs_sweep_ridge_bayes_cpp(arma::mat& Xc,
                                   arma::colvec& time, arma::colvec& delta,
                                   arma::imat& R,
                                   arma::colvec& beta,
                                   arma::vec& ch_times, arma::vec& ch_vals,
                                   double lambda_j,
                                   arma::ivec& miss_vars,
                                   arma::vec& cv_lam) {
  int n = Xc.n_rows, p = Xc.n_cols;

  for(int jj = 0; jj < (int)miss_vars.n_elem; jj++) {
    int j = miss_vars(jj);
    int pm1 = p - 1;

    std::vector<int> obs_idx, miss_idx;
    for(int i = 0; i < n; i++) {
      if(R(i, j) == 0) miss_idx.push_back(i);
      else obs_idx.push_back(i);
    }
    if(miss_idx.empty()) continue;
    int n_obs = (int)obs_idx.size();

    arma::mat Z_obs(n_obs, pm1);
    arma::colvec y_obs(n_obs);
    for(int ii = 0; ii < n_obs; ii++) {
      int i = obs_idx[ii];
      y_obs(ii) = Xc(i, j);
      int col = 0;
      for(int k = 0; k < p; k++) {
        if(k == j) continue;
        Z_obs(ii, col++) = Xc(i, k);
      }
    }

    double y_mean = arma::mean(y_obs);
    arma::colvec yc = y_obs - y_mean;
    arma::rowvec Z_mean = arma::mean(Z_obs, 0);
    arma::mat Zc = Z_obs.each_row() - Z_mean;

    double lambda_ridge = cv_lam(jj);

    arma::mat ZtZ = Zc.t() * Zc / (double)n_obs;
    ZtZ.diag() += lambda_ridge;
    arma::colvec Zty = Zc.t() * yc / (double)n_obs;
    arma::colvec gamma_hat;
    bool solved = arma::solve(gamma_hat, ZtZ, Zty);
    if(!solved) { ZtZ.diag() += lambda_ridge; arma::solve(gamma_hat, ZtZ, Zty); }

    double alpha_hat = y_mean - arma::dot(Z_mean, gamma_hat);

    arma::colvec resid_obs = y_obs - alpha_hat - Z_obs * gamma_hat;
    int df = std::max(n_obs - pm1 - 1, 1);
    double sigma2_hat = std::max(arma::dot(resid_obs, resid_obs) / (double)df, 0.001);

    // Draw sigma2 using R::rchisq
    double chi2_draw = R::rchisq((double)df);
    double sigma2_draw = (double)df * sigma2_hat / std::max(chi2_draw, 1e-10);
    sigma2_draw = std::max(sigma2_draw, 0.001);

    // Single Cholesky for both solve and draw
    arma::mat L_ZtZ;
    bool chol_ok = arma::chol(L_ZtZ, ZtZ, "lower");
    arma::colvec gamma_draw;
    if(chol_ok) {
      double scale = sqrt(sigma2_draw / (double)n_obs);
      arma::colvec z_vec(pm1);
      for(int dd = 0; dd < pm1; dd++) z_vec(dd) = R::rnorm(0.0, 1.0);
      gamma_draw = gamma_hat + scale * arma::solve(trimatu(L_ZtZ.t()), z_vec);
    } else {
      gamma_draw = gamma_hat;
    }

    double alpha_draw = y_mean - arma::dot(Z_mean, gamma_draw);
    double sigma_draw = sqrt(sigma2_draw);

    for(int ii = 0; ii < (int)miss_idx.size(); ii++) {
      int i = miss_idx[ii];

      arma::rowvec Xi_mj(pm1);
      int col = 0;
      for(int k = 0; k < p; k++) {
        if(k == j) continue;
        Xi_mj(col++) = Xc(i, k);
      }

      double mu = alpha_draw + arma::dot(Xi_mj, gamma_draw);
      double x_star = R::rnorm(mu, sigma_draw);
      x_star = std::max(-5.0, std::min(5.0, x_star));
      double x_curr = Xc(i, j);

      arma::rowvec Xi_star = Xc.row(i); Xi_star(j) = x_star;
      arma::rowvec Xi_curr = Xc.row(i); Xi_curr(j) = x_curr;

      double f_star = cox_lik_individual(time(i), delta(i), Xi_star, beta,
                                         ch_times, ch_vals);
      double f_curr = cox_lik_individual(time(i), delta(i), Xi_curr, beta,
                                         ch_times, ch_vals);

      double ratio = std::min(1.0, f_star / std::max(f_curr, 1e-30));
      if(R::runif(0, 1) < ratio) {
        Xc(i, j) = x_star;
      }
    }
  }
}

// ================================================================
//  Point-estimate ridge sweep (used by Standard IRO / Liang method)
//  Also optimized with R::rchisq and single Cholesky
// ================================================================
// [[Rcpp::export]]
void smc_fcs_sweep_cpp(arma::mat& Xc,
                       arma::colvec& time, arma::colvec& delta,
                       arma::imat& R,
                       arma::colvec& beta,
                       arma::vec& ch_times, arma::vec& ch_vals,
                       double lambda_j,
                       arma::ivec& miss_vars) {
  int n = Xc.n_rows, p = Xc.n_cols;
  for(int jj = 0; jj < (int)miss_vars.n_elem; jj++) {
    int j = miss_vars(jj);
    std::vector<int> miss_idx;
    for(int i = 0; i < n; i++) {
      if(R(i, j) == 0) miss_idx.push_back(i);
    }
    if(miss_idx.empty()) continue;
    arma::mat Z(n, p - 1); int col = 0;
    for(int k = 0; k < p; k++) { if(k == j) continue; Z.col(col++) = Xc.col(k); }
    arma::colvec y = Xc.col(j);
    int n_obs = 0; double y_sum = 0.0;
    arma::rowvec Z_sum(p-1, fill::zeros);
    for(int i = 0; i < n; i++) {
      if(R(i,j)==1) { n_obs++; y_sum += y(i); Z_sum += Z.row(i); }
    }
    double n_obs_d = std::max((double)n_obs, 1.0);
    double y_mean = y_sum / n_obs_d;
    arma::rowvec Z_mean = Z_sum / n_obs_d;
    arma::colvec yc = y - y_mean;
    arma::mat Zc = Z.each_row() - Z_mean;
    double lambda_ridge = std::max((double)(p-1) / (n_obs_d * n_obs_d), 1e-8);
    arma::mat ZtZ = Zc.t() * Zc / (double)n; ZtZ.diag() += lambda_ridge;
    arma::colvec Zty = Zc.t() * yc / (double)n;
    arma::colvec gamma_hat;
    bool solved = arma::solve(gamma_hat, ZtZ, Zty);
    if(!solved) { ZtZ.diag() += lambda_ridge; arma::solve(gamma_hat, ZtZ, Zty); }
    double alpha_hat = y_mean - arma::dot(Z_mean, gamma_hat);
    double ss_obs = 0.0;
    for(int i = 0; i < n; i++) {
      if(R(i,j)==1) { double r_i = y(i)-alpha_hat-arma::dot(Z.row(i),gamma_hat); ss_obs += r_i*r_i; }
    }
    double sigma = sqrt(std::max(ss_obs / n_obs_d, 0.01));
    for(int ii = 0; ii < (int)miss_idx.size(); ii++) {
      int i = miss_idx[ii];
      arma::rowvec Xi_mj(p-1); col = 0;
      for(int k = 0; k < p; k++) { if(k==j) continue; Xi_mj(col++) = Xc(i,k); }
      double mu = alpha_hat + arma::dot(Xi_mj, gamma_hat);
      double x_star = std::max(-5.0, std::min(5.0, R::rnorm(mu, sigma)));
      arma::rowvec Xi_star = Xc.row(i); Xi_star(j) = x_star;
      arma::rowvec Xi_curr = Xc.row(i);
      double f_star = cox_lik_individual(time(i),delta(i),Xi_star,beta,ch_times,ch_vals);
      double f_curr = cox_lik_individual(time(i),delta(i),Xi_curr,beta,ch_times,ch_vals);
      if(R::runif(0,1) < std::min(1.0, f_star/std::max(f_curr,1e-30)))
        Xc(i,j) = x_star;
    }
  }
}
