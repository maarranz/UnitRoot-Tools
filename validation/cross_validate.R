# Stage 2 validation only. No package code is installed or modified by this script.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop('Usage: Rscript cross_validate.R series.csv results.csv metadata.json [R-library]')
if (length(args) >= 4 && nzchar(args[4])) .libPaths(c(args[4], .libPaths()))
for (p in c('urca', 'bootUR', 'jsonlite')) if (!requireNamespace(p, quietly = TRUE)) stop(paste('Missing R package:', p))
if (as.character(utils::packageVersion('urca')) != '1.3.4' ||
    as.character(utils::packageVersion('bootUR')) != '1.0.5') stop('This audit requires urca 1.3-4 and bootUR 1.0.5')
options(digits = 17)
d <- read.csv(args[1], check.names = FALSE)
rows <- list()
add <- function(id, test, trend, lags, variant, statistic, pvalue = NA_real_, nobs = NA_real_, cv = rep(NA_real_, 3), estimate = NA_real_, reported_lags = lags) {
  rows[[length(rows) + 1L]] <<- data.frame(series = id, backend = 'R', test = test, trend = trend,
    lags = lags, variant = variant, statistic = as.numeric(statistic), pvalue = as.numeric(pvalue),
    nobs = nobs, cv1 = cv[1], cv5 = cv[2], cv10 = cv[3], estimate = as.numeric(estimate), reported_lags = reported_lags)
}
for (id in names(d)) {
  y <- d[[id]]; y <- y[!is.na(y)]
  for (tr in c('c', 'ct')) {
    model <- if (tr == 'c') 'constant' else 'trend'
    for (p in c(0L, 4L)) {
      f <- urca::ur.ers(y, type = 'DF-GLS', model = model, lag.max = p)
      add(id, 'DFGLS', tr, p, 'fixed', f@teststat, nobs = length(y) - 1L - p,
          cv = as.numeric(f@cval), estimate = f@testreg$coefficients[1, 1])
    }
    for (b in c(4L, 12L)) for (typ in c('tau', 'rho')) {
      f <- urca::ur.pp(y, type = if (typ == 'tau') 'Z-tau' else 'Z-alpha', model = model, use.lag = b)
      add(id, paste0('PP_', typ), tr, b, 'fixed', f@teststat, nobs = length(y) - 1L,
          cv = as.numeric(f@cval), estimate = f@testreg$coefficients[2, 1])
    }
  }
  for (tr in c('n', 'c', 'ct')) for (K in c(4L, 12L)) {
    det <- switch(tr, n = 'none', c = 'intercept', ct = 'trend')
    for (crit in c('MAIC', 'MBIC')) for (two in c(TRUE, FALSE)) {
      f <- bootUR::adf(y, deterministics = det, min_lag = 0L, max_lag = K,
                      criterion = crit, criterion_scale = FALSE, two_step = two)
      reported <- unname(f$details[['selected lags']])
      # bootUR 1.0.5 one-step reports parameter count minus one, including
      # deterministic coefficients; preserve raw output and label actual p.
      p <- reported - if (two) 0L else switch(tr, n = 0L, c = 1L, ct = 2L)
      variant <- paste0(if (two) 'two_step' else 'one_step', '_K', K, '_unscaled')
      add(id, crit, tr, p, variant, f$statistic, f$p.value, length(y) - 1L - p,
          estimate = f$estimate, reported_lags = reported)
    }
    # Sensitivity to bootUR's default volatility rescaling; not a matched ARCH specification.
    f <- bootUR::adf(y, deterministics = det, min_lag = 0L, max_lag = K,
                    criterion = 'MAIC', criterion_scale = TRUE, two_step = TRUE)
    p <- unname(f$details[['selected lags']])
    add(id, 'MAIC', tr, p, paste0('two_step_K', K, '_scaled'), f$statistic, f$p.value,
        length(y) - 1L - p, estimate = f$estimate)
  }
}
write.csv(do.call(rbind, rows), args[2], row.names = FALSE, na = '')
metadata <- list(R = R.version.string, platform = R.version$platform,
  packages = lapply(c('urca', 'bootUR', 'jsonlite', 'Rcpp', 'RcppArmadillo', 'RcppParallel', 'RcppThread'),
    function(p) list(name = p, version = as.character(utils::packageVersion(p)))),
  BLAS = unname(extSoftVersion()['BLAS']))
writeLines(jsonlite::toJSON(metadata, auto_unbox = TRUE, pretty = TRUE), args[3])
