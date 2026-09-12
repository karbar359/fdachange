#' Built-in RKHS kernels used throughout the dissertation.
#'
#' Each kernel is a function of two scalars.
kernel_gaussian <- function(x, y) exp(-(x - y)^2 / 2)
kernel_min      <- function(x, y) pmin(x, y)
kernel_2_max    <- function(x, y) 2 - pmax(x, y)

#' Evaluate a kernel on every pair (X_i, X_j) and return the Gram matrix.
#'
#' The built-in kernels (exp, pmin, pmax) are already vectorised, so they go
#' through `outer` directly: bit-identical to the old `Vectorize` wrapping
#' (verified for all three kernels) and ~90x faster, which is what makes the
#' n = 2000 null-cache regeneration tractable. A genuinely scalar-only kernel
#' still works via the fallback.
gram <- function(X, kernel) {
  out <- tryCatch(outer(X, X, FUN = kernel), error = function(e) NULL)
  if (is.matrix(out) && identical(dim(out), c(length(X), length(X)))) {
    return(out)
  }
  outer(X, X, FUN = Vectorize(kernel))
}

#' Resolve a kernel name ("gaussian" | "min" | "2_max") to its function.
kernel_fn <- function(name) {
  switch(name,
    gaussian = kernel_gaussian,
    min      = kernel_min,
    `2_max`  = kernel_2_max,
    stop("unknown kernel: ", name))
}
