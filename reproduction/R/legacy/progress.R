#' Format a duration in seconds as a short human-readable string.
#' Used by long-running simulation scripts to emit elapsed / ETA lines.
fmt_duration <- function(s) {
  if (!is.finite(s) || s < 0) return("--")
  if (s < 60)     return(sprintf("%.1fs", s))
  if (s < 3600)   return(sprintf("%dm%02ds", as.integer(s) %/% 60L, as.integer(s) %% 60L))
  if (s < 86400)  return(sprintf("%dh%02dm", as.integer(s) %/% 3600L, (as.integer(s) %% 3600L) %/% 60L))
  sprintf("%dd%02dh", as.integer(s) %/% 86400L, (as.integer(s) %% 86400L) %/% 3600L)
}

#' mclapply with in-cell progress reporting.
#'
#' Splits the item vector into batches of `batch_size` (default 4 * ncores)
#' and runs each batch through `parallel::mclapply` with
#' `mc.preschedule = TRUE` (so each batch is only `ncores` forks, not one
#' fork per item). Prints a single \r-overwriting line per batch with
#' count done / total / %, elapsed, and ETA.
#'
#' Matches the null-generator use case: dozens of cells, each running
#' 100-10000 MC replicates that should be visibly ticking while the
#' user waits.
mclapply_with_progress <- function(items, fn, ncores,
                                   prefix = "    ",
                                   batch_size = 4L * ncores) {
  n <- length(items)
  if (n == 0L) return(list())
  t0 <- Sys.time()
  out <- vector("list", n)
  done <- 0L
  while (done < n) {
    hi <- min(done + batch_size, n)
    chunk <- parallel::mclapply(
      items[(done + 1L):hi], fn,
      mc.cores = ncores, mc.set.seed = TRUE, mc.preschedule = TRUE)
    err <- which(vapply(chunk, inherits, logical(1), "try-error"))
    if (length(err)) {
      stop(sprintf("mclapply errored at item(s) %s",
                   paste(err + done, collapse = ", ")))
    }
    out[(done + 1L):hi] <- chunk
    done <- hi
    # Advance the parent's RNG state before the next batch so the child
    # fork substreams are distinct each round. Without this, every batch
    # re-uses the same N_CORES substreams and the output collapses onto
    # batch_size unique values (observed: 30000 draws -> 80 uniques when
    # generating the Wilcoxon gamma = 0 null before this line was added).
    invisible(runif(1))
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    eta <- if (done > 0L) elapsed / done * (n - done) else NA_real_
    cat(sprintf("\r%s%d/%d (%d%%)  elapsed=%s  eta=%s       ",
                prefix, done, n, as.integer(100 * done / n),
                fmt_duration(elapsed), fmt_duration(eta)))
    flush.console()
  }
  cat("\n"); flush.console()
  out
}
