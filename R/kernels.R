#' Built-in RKHS kernels
#'
#' The three reproducing kernels used by the dissertation's RKHS scan:
#' Gaussian `exp(-(x - y)^2 / 2)`, min `min(x, y)` (the Wiener kernel on
#' `[0, 1]`), and 2-max `2 - max(x, y)`. Each is a vectorised function of
#' two numeric arguments.
#'
#' @param x,y numeric vectors.
#' @return numeric vector of kernel evaluations.
#' @seealso [gram()], [rkhs_statistic()]
#' @examples
#' kernel_gaussian(0, c(-1, 0, 1))
#' kernel_min(0.3, 0.7)
#' kernel_2max(0.3, 0.7)
#' @export
kernel_gaussian <- function(x, y) exp(-(x - y)^2 / 2)

#' @rdname kernel_gaussian
#' @export
kernel_min <- function(x, y) pmin(x, y)

#' @rdname kernel_gaussian
#' @export
kernel_2max <- function(x, y) 2 - pmax(x, y)

#' Gram matrix of a sample under a kernel
#'
#' Evaluates `kernel` on every pair `(X_i, X_j)` of a numeric sample.
#' Vectorised kernels (all built-ins) go through [outer()] directly; a
#' genuinely scalar-only kernel still works via a `Vectorize` fallback.
#'
#' @param X numeric sample vector.
#' @param kernel a kernel function of two arguments, e.g.
#'   [kernel_gaussian()].
#' @return `length(X)` x `length(X)` numeric matrix.
#' @examples
#' gram(rnorm(5), kernel_gaussian)
#' @export
gram <- function(X, kernel) {
  out <- tryCatch(outer(X, X, FUN = kernel), error = function(e) NULL)
  if (is.matrix(out) && identical(dim(out), c(length(X), length(X)))) {
    return(out)
  }
  outer(X, X, FUN = Vectorize(kernel))
}

# Resolve a kernel given as a name ("gaussian" | "min" | "2max", with
# "2_max" and "2.max" accepted as legacy spellings) or as a function.
.resolve_kernel <- function(kernel) {
  if (is.function(kernel)) return(kernel)
  switch(as.character(kernel),
    gaussian = kernel_gaussian,
    min      = kernel_min,
    `2max`   = ,
    `2_max`  = ,
    `2.max`  = kernel_2max,
    stop("unknown kernel: ", kernel))
}

# Canonical short name of a kernel for cache keys ("gaussian", "min",
# "2max"). Functions other than the built-ins have no canonical name.
.kernel_name <- function(kernel) {
  if (is.character(kernel)) {
    k <- sub("^2[._]max$", "2max", kernel)
    if (k %in% c("gaussian", "min", "2max")) return(k)
    stop("unknown kernel name: ", kernel)
  }
  if (identical(kernel, kernel_gaussian)) return("gaussian")
  if (identical(kernel, kernel_min)) return("min")
  if (identical(kernel, kernel_2max)) return("2max")
  NA_character_
}
