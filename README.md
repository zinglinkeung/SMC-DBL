# SMC-DBL

Research code for Cox proportional hazards regression with missing covariates, combining substantive-model-compatible imputation, debiased lasso estimation, and Rubin pooling.

The project contains R simulation and analysis scripts backed by C++ routines for Cox likelihood calculations, Breslow baseline cumulative hazard estimation, and imputation sweeps. It also includes full-data, mean-imputation, complete-case, and iterative imputation comparison methods.


## Files

Keep these files together and run commands from their directory; paths in the scripts are relative to the current working directory.

| File | Purpose |
| --- | --- |
| `iro_smc_fcs_lib.cpp` | RcppArmadillo backend: Cox score, Hessian and score outer-product calculations; partial likelihood; Breslow estimator; ridge-based imputation sweeps; auxiliary Cox lasso and concordance routines. |
| `bench_single_iter.R` | Runs one simulation setting/replicate for full-data oracle, SMC-DBL, and the comparison labeled “Standard IRO (Liang).” Despite its filename, it runs full method workflows, not one update. |
| `bench_meanimp.R` | Runs column-mean imputation followed by debiased Cox lasso and appends the result to a replicate file. |
| `bench_complete.R` | Runs complete-case debiased Cox lasso and appends the result. Its introductory mean-imputation comment is stale; the implementation uses complete cases. |
| `collect_results.R` | Aggregates replicate files into bias, RMSE, coverage, variance calibration, and selection summaries. |
| `realdata_analysis.R` | Runs five methods on SNP and clinical predictors loaded from the external `df_all.csv`. |

## Dependencies

Install R and a working C++ toolchain compatible with your R installation. The analysis scripts compile the backend at startup using `Rcpp::sourceCpp()`; a precompiled binary is not supplied. Typical build environments are Rtools on Windows, the R-compatible compiler tools on macOS, or development compilers on Linux.

Install the packages explicitly loaded by the scripts:

```r
install.packages(
  c("Rcpp", "RcppArmadillo", "glmnet", "quadprog", "MASS", "survival", "smcfcs"),
  repos = "https://cloud.r-project.org"
)
```
