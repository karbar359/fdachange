# Null-draw storage: shipped caches in inst/extdata/<method>_null.rds, user
# extensions in null_cache_dir()/<method>_null.rds. Each file is a named
# list of numeric draw vectors keyed by .null_key(); user entries override
# shipped ones. Reads are memoised per session.

.null_env <- new.env(parent = emptyenv())

.fmt_num <- function(x) format(x, trim = TRUE, scientific = FALSE)

# Cell key shared by the accessors, save_null_draws() and the extdata
# build script (data-raw/build_null_caches.R).
.null_key <- function(method, gamma, beta = 0, p = 2, kernel = "gaussian",
                      baseline = "") {
  switch(method,
    wilcoxon = sprintf("g%s_b%s", .fmt_num(gamma), .fmt_num(beta)),
    lp       = sprintf("g%s_p%s", .fmt_num(gamma), .fmt_num(p)),
    rkhs     = sprintf("g%s_%s%s", .fmt_num(gamma), .kernel_name(kernel),
                       if (nzchar(baseline)) paste0("_", baseline) else ""),
    birbilas = sprintf("g%s", .fmt_num(gamma)),
    stop("unknown method: ", method))
}

.null_table <- function(method) {
  memo <- paste0("tbl_", method)
  if (!is.null(.null_env[[memo]])) return(.null_env[[memo]])
  tbl <- list()
  shipped <- system.file("extdata", paste0(method, "_null.rds"),
                         package = "fdachange")
  if (nzchar(shipped)) tbl <- readRDS(shipped)
  user_path <- file.path(null_cache_dir(), paste0(method, "_null.rds"))
  if (file.exists(user_path)) {
    user <- readRDS(user_path)
    tbl[names(user)] <- user
  }
  .null_env[[memo]] <- tbl
  tbl
}

.null_get <- function(method, key, simulate_hint) {
  tbl <- .null_table(method)
  draws <- tbl[[key]]
  if (is.null(draws)) {
    stop(sprintf(paste0(
      "No %s null draws for cell \"%s\". Available cells: %s.\n",
      "Generate the missing cell with %s and store it with save_null_draws()."),
      method, key,
      if (length(tbl)) paste(names(tbl), collapse = ", ") else "(none)",
      simulate_hint), call. = FALSE)
  }
  draws
}

#' Directory for user-generated null draws
#'
#' Null draws saved with [save_null_draws()] live in one `.rds` file per
#' method under this directory, and take precedence over the caches
#' shipped with the package. Defaults to
#' `tools::R_user_dir("fdachange", "cache")`; override for a session with
#' `options(fdachange.cache_dir = ...)`.
#'
#' @return path of the cache directory (not necessarily existing yet).
#' @examples
#' null_cache_dir()
#' @export
null_cache_dir <- function() {
  getOption("fdachange.cache_dir", tools::R_user_dir("fdachange", "cache"))
}

#' Store null draws for later lookup
#'
#' Saves a vector of Monte-Carlo null draws (typically from
#' [simulate_wilcoxon_null()] and friends) into [null_cache_dir()] under
#' the cell key that the corresponding `*_null_draws()` accessor uses, so
#' subsequent calls - including [epidemic_test()] - find it
#' automatically. An existing user entry for the same cell is
#' overwritten; shipped package caches are never modified.
#'
#' @param draws numeric vector of null draws.
#' @param method one of `"wilcoxon"`, `"lp"`, `"rkhs"`, `"birbilas"`.
#' @param gamma weight exponent of the cell.
#' @param beta Wilcoxon only: logarithmic weight exponent.
#' @param p Lp only: the Lp exponent.
#' @param kernel,baseline RKHS only: kernel name and H0-distribution tag,
#'   see [rkhs_null_draws()].
#' @return the cell key, invisibly.
#' @examples
#' \donttest{
#' draws <- simulate_wilcoxon_null(50, gamma = 0.1, grid = 200)
#' save_null_draws(draws, "wilcoxon", gamma = 0.1)
#' quantile(wilcoxon_null_draws(gamma = 0.1), 0.95)
#' }
#' @export
save_null_draws <- function(draws, method = c("wilcoxon", "lp", "rkhs",
                                              "birbilas"),
                            gamma, beta = 0, p = 2, kernel = "gaussian",
                            baseline = "") {
  method <- match.arg(method)
  stopifnot(is.numeric(draws), length(draws) > 0)
  key <- .null_key(method, gamma = gamma, beta = beta, p = p,
                   kernel = kernel, baseline = baseline)
  dir <- null_cache_dir()
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(dir, paste0(method, "_null.rds"))
  user <- if (file.exists(path)) readRDS(path) else list()
  user[[key]] <- as.numeric(draws)
  saveRDS(user, path, compress = "xz")
  rm(list = ls(envir = .null_env), envir = .null_env)  # drop memoised tables
  invisible(key)
}
