#!/usr/bin/env Rscript
# Fetch the raw LITGRID balancing-energy-price CSV exports used by the
# application scripts, into <repo>/data-raw/litgrid-raw/ (gitignored,
# ~119 MB unpacked).
#
# Primary source: the GitHub release asset `litgrid-raw-v1.tar.gz` of this
# repository - the exact 20 exports the dissertation analysed, frozen. The
# data originates from the public Baltic transparency dashboard
# (https://baltic.transparency-dashboard.eu/), but the dashboard serves
# live exports whose history can be revised, so reproduction pins the
# frozen copies and verifies every file against the committed SHA256SUMS.
#
# Usage:  Rscript reproduction/applications/litgrid/download_litgrid.R
# Env:    LITGRID_RAW_DIR overrides the target directory.
#         GITHUB_REPO overrides the release source (default
#         karbar359/fdachange; a private repo needs a GITHUB_TOKEN with
#         repo scope, or use `gh release download` manually, see below).

RELEASE_TAG <- "litgrid-raw-v1"
ASSET       <- "litgrid-raw-v1.tar.gz"
REPO        <- Sys.getenv("GITHUB_REPO", "karbar359/fdachange")

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- sub("--file=", "", args[grepl("--file=", args)])
  if (length(hit)) dirname(normalizePath(hit[1])) else getwd()
})
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
target <- Sys.getenv("LITGRID_RAW_DIR",
                     file.path(repo_root, "data-raw", "litgrid-raw"))
dir.create(target, showWarnings = FALSE, recursive = TRUE)

sums_file <- file.path(script_dir, "SHA256SUMS")
stopifnot(file.exists(sums_file))
sums <- read.table(sums_file, col.names = c("sha256", "file"),
                   stringsAsFactors = FALSE)

verify <- function() {
  ok <- TRUE
  for (i in seq_len(nrow(sums))) {
    f <- file.path(target, sums$file[i])
    if (!file.exists(f)) { cat("MISSING ", sums$file[i], "\n"); ok <- FALSE; next }
    dig <- system2("sha256sum", shQuote(f), stdout = TRUE)
    dig <- strsplit(dig, "\\s+")[[1]][1]
    if (!identical(dig, sums$sha256[i])) {
      cat("CHECKSUM MISMATCH ", sums$file[i], "\n"); ok <- FALSE
    }
  }
  ok
}

if (verify()) {
  cat("All", nrow(sums), "files already present and verified in", target, "\n")
  quit(status = 0)
}

tarball <- file.path(tempdir(), ASSET)
gh <- Sys.which("gh")
if (nzchar(gh)) {
  cat("Downloading release asset via gh ...\n")
  status <- system2(gh, c("release", "download", RELEASE_TAG,
                          "--repo", REPO, "--pattern", ASSET,
                          "--dir", tempdir(), "--clobber"))
  if (status != 0) stop("gh release download failed")
} else {
  url <- sprintf("https://github.com/%s/releases/download/%s/%s",
                 REPO, RELEASE_TAG, ASSET)
  cat("Downloading", url, "\n(private repos need gh or a token; see header)\n")
  utils::download.file(url, tarball, mode = "wb", quiet = FALSE)
}
cat("Unpacking into", target, "\n")
utils::untar(tarball, exdir = target)

if (!verify()) stop("Checksum verification FAILED after download.")
cat("OK:", nrow(sums), "files verified against SHA256SUMS.\n")
