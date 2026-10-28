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
#
# RNG: every item draws from its own L'Ecuyer-CMRG substream (see
# .item_streams), so the draws are independent, identical for any `cores`,
# and reproducible under set.seed() whatever generator the caller uses.
# mclapply's own mc.set.seed is not used: it resets its stream base to the
# parent's .Random.seed on every call and never advances the parent, so a
# batched loop re-forks the same substreams each round (and advancing the
# parent by one uniform between batches only shifts every second batch by
# one normal draw, which is what produced the dependent null caches).
.apply_progress <- function(items, fn, cores = 1L, progress = FALSE,
                            prefix = "  ") {
  n <- length(items)
  if (n == 0L) return(list())
  cores <- max(1L, as.integer(cores))
  batch_size <- if (cores > 1L) 4L * cores else max(1L, n %/% 20L)
  rng <- .item_streams(n)
  on.exit(assign(".Random.seed", rng$after, envir = globalenv()), add = TRUE)
  run <- function(i) {
    assign(".Random.seed", rng$streams[[i]], envir = globalenv())
    fn(items[[i]])
  }
  t0 <- Sys.time()
  out <- vector("list", n)
  done <- 0L
  while (done < n) {
    hi <- min(done + batch_size, n)
    idx <- (done + 1L):hi
    chunk <- if (cores > 1L) {
      parallel::mclapply(idx, run, mc.cores = cores, mc.set.seed = FALSE,
                         mc.preschedule = TRUE)
    } else {
      lapply(idx, run)
    }
    err <- which(vapply(chunk, inherits, logical(1), "try-error"))
    if (length(err)) {
      stop(sprintf("parallel worker errored at item(s) %s",
                   paste(err + done, collapse = ", ")))
    }
    out[idx] <- chunk
    done <- hi
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

# One L'Ecuyer-CMRG seed per item, chained with parallel::nextRNGStream.
# An L'Ecuyer caller seeds the chain from its current state and is left on
# the stream after the last item. Any other caller spends one draw of its
# own generator on the chain's seed and gets its state (and kind: the first
# element of .Random.seed encodes it) back afterwards. Returns the streams
# and the .Random.seed to restore once the items have run.
.item_streams <- function(n) {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1)
  }
  lecuyer <- RNGkind()[1L] == "L'Ecuyer-CMRG"
  if (lecuyer) {
    s <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    seed <- sample.int(.Machine$integer.max, 1L)
    saved <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    s <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    assign(".Random.seed", saved, envir = globalenv())
  }
  streams <- vector("list", n)
  for (i in seq_len(n)) {
    streams[[i]] <- s
    s <- parallel::nextRNGStream(s)
  }
  after <- if (lecuyer) s else saved
  list(streams = streams, after = after)
}
