#!/usr/bin/env Rscript
# Descriptive tail/spike statistics for the Litgrid balancing-price series,
# backing the statement "The prices are heavy tailed and spiky."
#
# Computes, side by side for two versions of the hourly series:
#   analysed  - the exact 344 x 168 weekly matrix the tests and figures use
#               (load_litgrid_balancing_prices + the run.R trim), flattened.
#   dedup_raw - the full raw hourly series with each hour label counted ONCE.
#               The loader sums rows sharing an hour label, and the raw
#               exports overlap (the 2018 file runs to 2019-10-01, the 2019
#               file starts 2019-01-01), so ~6.5k hours in Jan-Sep 2019 are
#               doubled in the analysed matrix (DST fall-back hours too).
#               This column shows the conclusions do not hinge on that.
#
# Statistics: moments (skewness, excess kurtosis), quantiles, Gaussian
# exceedance benchmark (share beyond 4/6 sd vs the normal expectation),
# Hill tail index on |x| at two k, spike-episode run lengths above the
# 99.5% |x| quantile, and first-difference kurtosis / largest hourly jump.
#
# Results CSV goes to data/ (committed); runtime is a few seconds.
#
# Usage:  Rscript reproduction/litgrid_tail_stats.R
# Env:    LITGRID_RAW_DIR overrides the raw-data lookup (default:
#         this script until it finds a .raw sibling).

# ---- Inputs (edit these) --------------------------------------------------
ROW_DROP    <- 204L                 # outlier week (matches applications/litgrid/run.R)
ROW_KEEP_N  <- 344L
Q_PROBS     <- c(0.001, 0.01, 0.25, 0.75, 0.99, 0.999)
HILL_FRACS  <- c(0.025, 0.05)       # top-|x| fractions k/n for the Hill estimator
SPIKE_Q     <- 0.995                # |x| quantile defining a spike hour
SIGMA_MULTS <- c(4, 6)              # Gaussian-benchmark exceedance thresholds
OUT_CSV     <- "litgrid_tail_stats.csv"   # written under data/
# ---------------------------------------------------------------------------

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  key <- "--file="
  hit <- sub(key, "", args[grepl(key, args)])
  if (length(hit)) return(normalizePath(hit[1]))
  NULL
})()
root <- if (is.null(script_path)) getwd() else dirname(dirname(script_path))

source(file.path(root, "applications", "litgrid", "data.R"))

litgrid_dir <- Sys.getenv("LITGRID_RAW_DIR", unset = "")
if (!nzchar(litgrid_dir)) {
  litgrid_dir <- normalizePath(file.path(root, "..", "data-raw", "litgrid-raw"),
                               mustWork = FALSE)
}
if (!dir.exists(litgrid_dir)) {
  stop("Litgrid raw CSVs not found at ", litgrid_dir,
       "; run Rscript applications/litgrid/download_litgrid.R or set LITGRID_RAW_DIR.")
}
litgrid_dir <- litgrid_dir

# ---- Series 1: the analysed weekly matrix, flattened -----------------------
f <- load_litgrid_balancing_prices(litgrid_dir)
kept <- setdiff(seq_len(nrow(f)), ROW_DROP)
f <- f[kept, , drop = FALSE][seq_len(ROW_KEEP_N), , drop = FALSE]
x_analysed  <- as.numeric(t(f))                 # chronological: wk1 h1..168, wk2, ...
dx_analysed <- as.numeric(apply(f, 1, diff))    # within-week hour-to-hour jumps
cat(sprintf("analysed : %d weeks x %d hours = %d values, weeks %s..%s\n",
            nrow(f), ncol(f), length(x_analysed),
            rownames(f)[1], rownames(f)[nrow(f)]))

# ---- Series 2: raw hourly series, each hour label counted once -------------
files <- list.files(litgrid_dir,
                    pattern = "^balancing-energy-prices.*\\.csv$",
                    full.names = TRUE)
raw <- readr::read_delim(files, delim = ";", escape_double = FALSE,
                         trim_ws = TRUE, id = "id", show_col_types = FALSE)
d <- subset(raw, area == "Lithuania")
d$value    <- as.numeric(gsub(",", ".", d$Upward)) -
              as.numeric(gsub(",", ".", d$Downward))
d$datetime <- lubridate::floor_date(lubridate::ymd_hm(d$datetime_to), "hour")
d <- d[order(d$datetime), ]
n_dup <- sum(duplicated(d$datetime))
d <- d[!duplicated(d$datetime), ]
x_dedup  <- d$value
dx_dedup <- diff(x_dedup)
cat(sprintf("dedup_raw: %d hourly values (%d duplicated hour labels removed), %s..%s\n",
            length(x_dedup), n_dup,
            format(min(d$datetime)), format(max(d$datetime))))

# ---- Statistics -------------------------------------------------------------
tail_stats <- function(x, dx) {
  n  <- length(x)
  m  <- mean(x); s <- sd(x)
  m2 <- mean((x - m)^2)
  qs <- quantile(x, probs = Q_PROBS, names = FALSE)
  xs <- sort(abs(x), decreasing = TRUE)
  hill <- vapply(HILL_FRACS, function(fr) {
    k <- ceiling(fr * n)
    k / sum(log(xs[1:k] / xs[k + 1]))
  }, numeric(1))
  exceed <- vapply(SIGMA_MULTS, function(mu) sum(abs(x - m) > mu * s), numeric(1))
  expect <- vapply(SIGMA_MULTS, function(mu) n * 2 * pnorm(-mu), numeric(1))
  thr  <- quantile(abs(x), SPIKE_Q, names = FALSE)
  r    <- rle(abs(x) > thr)
  runs <- r$lengths[r$values]
  dm <- mean(dx); dm2 <- mean((dx - dm)^2)
  c(n_hours              = n,
    zero_share           = mean(x == 0),
    mean                 = m,
    sd                   = s,
    median               = median(x),
    iqr                  = IQR(x),
    skewness             = mean((x - m)^3) / m2^1.5,
    excess_kurtosis      = mean((x - m)^4) / m2^2 - 3,
    q0.1pct              = qs[1],
    q1pct                = qs[2],
    q25pct               = qs[3],
    q75pct               = qs[4],
    q99pct               = qs[5],
    q99.9pct             = qs[6],
    min                  = min(x),
    max                  = max(x),
    over4sd_hours        = exceed[1],
    over4sd_gauss_expect = expect[1],
    over6sd_hours        = exceed[2],
    over6sd_gauss_expect = expect[2],
    hill_alpha_k2.5pct   = hill[1],
    hill_alpha_k5pct     = hill[2],
    spike_thr_abs        = thr,
    spike_episodes       = length(runs),
    spike_len_median_h   = median(runs),
    spike_len_max_h      = max(runs),
    diff_excess_kurtosis = mean((dx - dm)^4) / dm2^2 - 3,
    max_abs_hourly_jump  = max(abs(dx)))
}

s_a <- tail_stats(x_analysed, dx_analysed)
s_d <- tail_stats(x_dedup,   dx_dedup)

out <- data.frame(statistic = names(s_a),
                  analysed  = as.character(signif(unname(s_a), 6)),
                  dedup_raw = as.character(signif(unname(s_d), 6)))

cat("\n")
print(out, row.names = FALSE, right = FALSE)

out_path <- file.path(root, "data", OUT_CSV)
write.csv(out, out_path, row.names = FALSE, quote = FALSE)
cat(sprintf("\nWrote %s\n", out_path))
