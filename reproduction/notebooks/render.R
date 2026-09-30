#!/usr/bin/env Rscript
# Render every reproduction notebook to its committed github_document .md
# (plus the figure PNGs under <name>_files/). Run from anywhere:
#
#   Rscript reproduction/notebooks/render.R
#
# The LITGRID notebook's raw-data sections evaluate only when the raw
# download is present (applications/litgrid/download_litgrid.R). The
# committed .md was rendered WITH the raw data; rendering without it
# still succeeds, with those sections skipped.

here <- dirname(normalizePath(sub("--file=", "",
  grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))

if (!nzchar(Sys.which("pandoc"))) {   # local installs live in ~/.local/bin
  Sys.setenv(PATH = paste(file.path(Sys.getenv("HOME"), ".local", "bin"),
                          Sys.getenv("PATH"), sep = .Platform$path.sep))
}

for (rmd in list.files(here, pattern = "[.]Rmd$", full.names = TRUE)) {
  cat("== rendering", basename(rmd), "\n")
  rmarkdown::render(rmd, output_format = "github_document",
                    envir = new.env(), quiet = TRUE)
}
cat("done\n")
