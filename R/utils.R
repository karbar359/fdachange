# Format a duration in seconds as a short human-readable string.
.fmt_duration <- function(s) {
  if (!is.finite(s) || s < 0) return("--")
  if (s < 60)     return(sprintf("%.1fs", s))
  if (s < 3600)   return(sprintf("%dm%02ds", as.integer(s) %/% 60L, as.integer(s) %% 60L))
  if (s < 86400)  return(sprintf("%dh%02dm", as.integer(s) %/% 3600L, (as.integer(s) %% 3600L) %/% 60L))
  sprintf("%dd%02dh", as.integer(s) %/% 86400L, (as.integer(s) %% 86400L) %/% 3600L)
}

# Apply fn over items, optionally in parallel (parallel::mclapply, batched
# so each batch is only `cores` forks) and optionally with a
# \r-overwriting progress line (count done / total / %, elapsed, ETA).
# Between batches the parent's RNG state advances so child fork substreams
# are distinct each round; without this every batch reuses the same
# substreams and the output collapses onto batch_size unique values.
.apply_progress <- function(items, fn, cores = 1L, progress = FALSE,
                            prefix = "  ") {
  n <- length(items)
  if (n == 0L) return(list())
  cores <- max(1L, as.integer(cores))
  batch_size <- if (cores > 1L) 4L * cores else max(1L, n %/% 20L)
  t0 <- Sys.time()
  out <- vector("list", n)
  done <- 0L
  while (done < n) {
    hi <- min(done + batch_size, n)
    chunk <- if (cores > 1L) {
      parallel::mclapply(items[(done + 1L):hi], fn, mc.cores = cores,
                         mc.set.seed = TRUE, mc.preschedule = TRUE)
    } else {
      lapply(items[(done + 1L):hi], fn)
    }
    err <- which(vapply(chunk, inherits, logical(1), "try-error"))
    if (length(err)) {
      stop(sprintf("parallel worker errored at item(s) %s",
                   paste(err + done, collapse = ", ")))
    }
    out[(done + 1L):hi] <- chunk
    done <- hi
    if (cores > 1L) invisible(runif(1))
    if (progress) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      eta <- if (done > 0L) elapsed / done * (n - done) else NA_real_
      cat(sprintf("\r%s%d/%d (%d%%)  elapsed=%s  eta=%s       ",
                  prefix, done, n, as.integer(100 * done / n),
                  .fmt_duration(elapsed), .fmt_duration(eta)))
      flush.console()
    }
  }
  if (progress) { cat("\n"); flush.console() }
  out
}
