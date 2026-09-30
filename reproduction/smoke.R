#!/usr/bin/env Rscript
# Cross-product smoke test for every (method, scenario) cell, ported onto
# the fdachange package (scenarios / projections / split aggregation) with
# the LEGACY statistics and null loaders from R/legacy/, so the table
# reproduces the original pipeline's smoke numbers:
#   methods   : rkhs | lp | wilcoxon
#   scenarios : S1 (normal) .. S5 (Pareto shape) -- scalar one-dim series
#               F1..F3 functional (paper 1 S1-S3) -- matrix(l+1, n), projected
#                                                    via c_max_var / ito / l2
#                                                    before the scalar test
# Prints a table of power at a = a_LO and a = a_HI (scenario-specific) over
# NREP replicates. Intent is to catch regressions end-to-end, not to measure
# calibration precisely.
#
# RKHS is pinned to the (kernel, baseline) pairs paper 2 studied, because the
# null cache only covers those. Wilcoxon runs at n = 200 because its
# asymptotic Hoelder-type null gives ~0 power at n = 50 regardless of a.
# Functional scenarios use the norm-based Wilcoxon kernel (h_0) because the
# changes are scale/tail shifts that the default median-based kernel (h_w)
# is invariant to; each functional scenario is run against all three paper-1
# projections.
#
# Usage: Rscript reproduction/smoke.R
# REPRO_CACHE_DIR may be set to override the cache directory lookup.

# ---- Inputs (edit these) --------------------------------------------------
GAMMA_GRID <- c(0, 0.25, 0.45)
NREP       <- 20
ALPHA      <- 0.05
EC         <- c(0.4, 0.7)
N_LP       <- 50
N_RKHS     <- 50
N_WILCOXON <- 200
N_F1       <- 200                  # functional cells at m = 1
N_F_M2     <- 300                  # functional cells at m = 2
N_STAB     <- 80                   # stability cell: small fixed-X sample
N_PROJ_STAB <- 20                  # stability cell: Ito draws per fixed X
NCORES     <- max(1L, parallel::detectCores() - 2L)
A_HI <- list(S1 = 2, S2 = 4, S3 = 2, S4 = 8, S5 = 1, F1 = 10, F2 = 5, F3 = 10, F4 = 2,
             F5 = 10)
A_LO <- list(S1 = 0, S2 = 0, S3 = 0, S4 = 0, S5 = 0, F1 = 2.5, F2 = 1.5, F3 = 2.5, F4 = 0,
             F5 = 0)
RKHS_PAIR <- list(
  S1 = list(kernel = "gaussian", baseline = ""),
  S2 = list(kernel = "min",      baseline = "rnorm.sq"),
  S4 = list(kernel = "2_max",    baseline = "rbeta.1.9"))
# ---------------------------------------------------------------------------

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  key <- "--file="
  hit <- sub(key, "", args[grepl(key, args)])
  if (length(hit)) return(normalizePath(hit[1]))
  ofile <- sys.frame(1)$ofile
  if (!is.null(ofile)) return(normalizePath(ofile))
  NULL
})()
root <- if (is.null(script_path)) getwd() else dirname(script_path)

suppressPackageStartupMessages(library(fdachange))
source(file.path(root, "R", "compat.R"))
source(file.path(root, "R", "legacy_env.R"))
legacy <- legacy_env(root)

project_identity <- projection_fn("none")
PC_D_N    <- 3L
FUNC_PROJ <- c("c_max_var", "ito", "l2", "wiener", "pc1", "pc_d")
resolve_projection <- function(name) {
  if (name == "pc_d") function(X) project_pc_d(X, d = PC_D_N)
  else               projection_fn(name)
}

# Tee all subsequent cat() output to results/smoke.txt so the table stays
# inspectable after the run (and between runs).
LOG_FILE <- file.path(root, "results", "smoke.txt")
dir.create(dirname(LOG_FILE), showWarnings = FALSE, recursive = TRUE)
sink(LOG_FILE, split = TRUE)
cat(sprintf("# %s  (REPRO_CACHE_DIR=%s, NCORES=%d)\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            Sys.getenv("REPRO_CACHE_DIR", unset = "<auto>"),
            NCORES))
# L'Ecuyer-CMRG gives reproducible parallel RNG streams across mclapply forks.
RNGkind("L'Ecuyer-CMRG")

cells <- list()
for (sc in names(A_HI)) {
  is_functional <- startsWith(sc, "F")

  if (!is_functional) {
    # Lp
    cells[[length(cells) + 1L]] <- list(
      method = "lp", sc = sc, N = N_LP, project = project_identity,
      stat_args = list(p = 2, cdf = get(paste0("cdf_", sc))),
      null_args = list(p = 2),
      note = "p=2")
    # RKHS
    if (sc %in% names(RKHS_PAIR)) {
      rp <- RKHS_PAIR[[sc]]
      cells[[length(cells) + 1L]] <- list(
        method = "rkhs", sc = sc, N = N_RKHS, project = project_identity,
        stat_args = list(kernel = kernel_fn(rp$kernel)),
        null_args = list(kernel_name = rp$kernel, baseline = rp$baseline),
        note = sprintf("%s/%s", rp$kernel,
                       if (nzchar(rp$baseline)) rp$baseline else "rnorm"))
    } else {
      cells[[length(cells) + 1L]] <- list(
        method = "rkhs", sc = sc, skip = "no null cache")
    }
    # Wilcoxon
    cells[[length(cells) + 1L]] <- list(
      method = "wilcoxon", sc = sc, N = N_WILCOXON, project = project_identity,
      stat_args = list(beta = 0, kernel = "h_w"),
      null_args = list(beta = 0),
      note = "h_w")
  } else {
    # Functional: Lp(ecdf) + Wilcoxon(h_0) at each projection.
    for (m_splits in c(1L, 2L)) {
      n_f <- if (m_splits == 1L) N_F1 else N_F_M2
      tag_m <- if (m_splits == 1L) "" else sprintf(" m=%d", m_splits)
      for (pname in FUNC_PROJ) {
        proj <- resolve_projection(pname)
        mult <- projection_multiplicity(pname, PC_D_N)
        cells[[length(cells) + 1L]] <- list(
          method = "lp", sc = sc, N = n_f, project = proj, M = m_splits,
          mult = mult,
          stat_args = list(p = 2, cdf = get(paste0("cdf_", sc))),
          null_args = list(p = 2),
          note = sprintf("p=2 x %s (ecdf)%s", pname, tag_m))
        cells[[length(cells) + 1L]] <- list(
          method = "wilcoxon", sc = sc, N = n_f, project = proj, M = m_splits,
          mult = mult,
          stat_args = list(beta = 0, kernel = "h_0"),
          null_args = list(beta = 0),
          note = sprintf("h_0 x %s%s", pname, tag_m))
      }
    }
    if (sc == "F4") {
      # F4 RKHS: hand-picked (projection, kernel) pairs with finite-null
      # cache (scenario_id = "F4" -> load_rkhs_F4_null).
      f4_rkhs_pairs <- list(
        c(projection = "c_max_var", kernel = "gaussian"),
        c(projection = "ito",       kernel = "gaussian"),
        c(projection = "l2",        kernel = "min"))
      for (pair in f4_rkhs_pairs) {
        pname <- pair[["projection"]]; kname <- pair[["kernel"]]
        cells[[length(cells) + 1L]] <- list(
          method = "rkhs", sc = sc, N = N_F1,
          project = resolve_projection(pname),
          stat_args = list(kernel = kernel_fn(kname)),
          null_args = list(kernel_name = kname,
                           scenario_id = "F4", projection = pname),
          note = sprintf("%s x %s", kname, pname))
      }
    } else {
      cells[[length(cells) + 1L]] <- list(
        method = "rkhs", sc = sc, skip = "no scenario null cache")
    }
    # Birbilas's E_{n,beta,gamma} consumes the full (l+1) x n matrix and is
    # only meaningful on FAR(1) scenarios. Wire on F4 only.
    if (sc == "F4") {
      cells[[length(cells) + 1L]] <- list(
        method = "far1", sc = sc, N = N_F1, project = project_identity,
        stat_args = list(beta_t = function(t) 0.7 * sqrt(t)),
        null_args = list(),
        matrix_in = TRUE,
        note = "E_n,beta,gamma (b=0.7sqrt(t))")
    }
  }
}

# Stability of random projections: fix one sample X and repeat the
# Ito-projection Wilcoxon scan N_PROJ_STAB times.
cells[[length(cells) + 1L]] <- list(
  method = "stability", sc = "F1", N = N_STAB,
  N_PROJ = N_PROJ_STAB,
  stat_args = list(beta = 0, kernel = "h_0"),
  note = sprintf("ito x N_proj=%d", N_PROJ_STAB))

cat(sprintf("%-9s %-3s %4s %-28s %5s %7s %7s\n",
            "method", "sc", "n", "note", "gamma", "p(a=LO)", "p(a=HI)"))
cat(strrep("-", 72), "\n", sep = "")

# Pre-warm raw-null CSV cache in the parent so every mclapply fork inherits
# the reads rather than hitting disk per worker.
invisible(tryCatch(legacy$raw_read("T2_gamma_beta.csv"),     error = function(e) NULL))
invisible(tryCatch(legacy$raw_read("dist_p2_ws_n200_t.csv"), error = function(e) NULL))
invisible(tryCatch(legacy$raw_read("SD_h_n2000_dist_t.csv"), error = function(e) NULL))
invisible(tryCatch(legacy$raw_read("E_q_gamma.csv"),         error = function(e) NULL))

run_cell <- function(cell) {
  if (!is.null(cell$skip)) {
    return(list(skip = TRUE, meta = cell))
  }
  if (cell$method == "stability") return(run_stability_cell(cell))
  m_splits <- if (is.null(cell$M))    1L else cell$M
  mult     <- if (is.null(cell$mult)) 1L else cell$mult
  stat_fn  <- get(paste0(cell$method, "_statistic"), envir = legacy)
  null_fn  <- get(paste0(cell$method, "_null_distribution"), envir = legacy)
  gen      <- get(paste0("generate_", cell$sc))
  pw_grid <- vector("list", length(GAMMA_GRID))
  for (gi in seq_along(GAMMA_GRID)) {
    g    <- GAMMA_GRID[gi]
    null <- do.call(null_fn, c(list(cell$N %/% m_splits, g), cell$null_args))
    cv   <- critical_value(null, ALPHA, splits = m_splits * mult)
    run  <- make_run(cell$project, stat_fn, g, cell$stat_args,
                     matrix_in = isTRUE(cell$matrix_in))
    pw <- numeric(2)
    for (ia in 1:2) {
      a <- if (ia == 1) A_LO[[cell$sc]] else A_HI[[cell$sc]]
      set.seed(42)
      stats <- numeric(NREP)
      for (i in seq_len(NREP)) {
        X <- gen(cell$N, a, EC)
        stats[i] <- scan_over_splits(X, m_splits, run)
      }
      pw[ia] <- mean(stats > cv)
    }
    pw_grid[[gi]] <- pw
  }
  list(skip = FALSE, meta = cell, pw_grid = pw_grid)
}

# Stability cell: fix X, draw N_PROJ Ito paths, compute T_n on each, return
# R_hat at a=LO and a=HI per gamma.
run_stability_cell <- function(cell) {
  gen <- get(paste0("generate_", cell$sc))
  pw_grid <- vector("list", length(GAMMA_GRID))
  for (gi in seq_along(GAMMA_GRID)) {
    g  <- GAMMA_GRID[gi]
    cv <- legacy$wilcoxon_critical_value(
      cell$N, g,
      beta  = cell$stat_args$beta,
      alpha = ALPHA, m.d = 1L)
    pw <- numeric(2)
    for (ia in 1:2) {
      a <- if (ia == 1) A_LO[[cell$sc]] else A_HI[[cell$sc]]
      set.seed(42)
      X  <- gen(cell$N, a, EC)              # FIXED for this (a, gamma)
      Ts <- numeric(cell$N_PROJ)
      for (i in seq_len(cell$N_PROJ)) {
        Xp    <- project_ito(X)
        Ts[i] <- do.call(legacy$wilcoxon_statistic,
                         c(list(Xp, g), cell$stat_args))
      }
      pw[ia] <- mean(Ts > cv)
    }
    pw_grid[[gi]] <- pw
  }
  list(skip = FALSE, meta = cell, pw_grid = pw_grid)
}

# mc.preschedule = FALSE gives each cell its own fork slot (work stealing).
results <- parallel::mclapply(cells, run_cell,
                              mc.cores = NCORES, mc.set.seed = TRUE,
                              mc.preschedule = FALSE)
if (any(vapply(results, inherits, logical(1), "try-error"))) {
  bad <- which(vapply(results, inherits, logical(1), "try-error"))
  stop(sprintf("mclapply errored in %d cell(s): %s",
               length(bad), paste(bad, collapse = ", ")))
}

for (res in results) {
  c_ <- res$meta
  if (isTRUE(res$skip)) {
    cat(sprintf("%-9s %-3s %4s %-28s %5s %7s %7s\n",
                c_$method, c_$sc, "-", c_$skip, "-", "-", "-"))
  } else {
    for (gi in seq_along(GAMMA_GRID)) {
      g  <- GAMMA_GRID[gi]
      pw <- res$pw_grid[[gi]]
      cat(sprintf("%-9s %-3s %4d %-28s %5.2f %7.2f %7.2f\n",
                  c_$method, c_$sc, c_$N, c_$note, g, pw[1], pw[2]))
    }
  }
}

sink()
message(sprintf("wrote %s", LOG_FILE))
