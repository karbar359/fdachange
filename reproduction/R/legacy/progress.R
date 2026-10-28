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
#' RNG: every item draws from its own L'Ecuyer-CMRG substream, chained with
#' `parallel::nextRNGStream` from the parent's current `.Random.seed` (the
#' scripts set `RNGkind("L'Ecuyer-CMRG")` and a seed first), and the parent
#' is left on the stream after the last item. Draws are therefore
#' independent and identical for any `ncores` and `batch_size`.
#'
#' `mc.set.seed = TRUE` is deliberately not used: mclapply resets its stream
#' base to the parent's `.Random.seed` on every call and never advances the
#' parent, so a batched loop re-forks the same substreams each round. The
#' earlier `runif(1)` between batches moved the parent by one uniform, half
#' a normal draw, so every second batch repeated an earlier batch's streams
#' shifted by one normal: 30,000 "draws" of the gamma = 0 Wilcoxon null held
#' about 160 distinct bridges (acf 0.994 at lag 160), and the RKHS and far1
#' caches of 1,000 held about 100.
mclapply_with_progress <- function(items, fn, ncores,
                                   prefix = "    ",
                                   batch_size = 4L * ncores) {
  n <- length(items)
  if (n == 0L) return(list())
  if (RNGkind()[1L] != "L'Ecuyer-CMRG") {
    stop("mclapply_with_progress needs RNGkind(\"L'Ecuyer-CMRG\") and a seed")
  }
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1)
  }
  s <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  streams <- vector("list", n)
  for (i in seq_len(n)) {
    streams[[i]] <- s
    s <- parallel::nextRNGStream(s)
  }
  assign(".Random.seed", s, envir = globalenv())
  run <- function(i) {
    assign(".Random.seed", streams[[i]], envir = globalenv())
    fn(items[[i]])
  }
  t0 <- Sys.time()
  out <- vector("list", n)
  done <- 0L
  while (done < n) {
    hi <- min(done + batch_size, n)
    idx <- (done + 1L):hi
    chunk <- parallel::mclapply(idx, run, mc.cores = ncores,
                                mc.set.seed = FALSE, mc.preschedule = TRUE)
    err <- which(vapply(chunk, inherits, logical(1), "try-error"))
    if (length(err)) {
      stop(sprintf("mclapply errored at item(s) %s",
                   paste(err + done, collapse = ", ")))
    }
    out[idx] <- chunk
    done <- hi
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
