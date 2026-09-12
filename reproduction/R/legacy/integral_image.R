#' Summed-area table over a 2D numeric matrix.
#'
#' @param img numeric matrix.
#' @return list with class "NumMatrix": `sat` (integral image, (r+1) x (c+1))
#'   and `empty` flag.
NumMatrix <- function(img) {
  if (!is.matrix(img) || nrow(img) == 0 || ncol(img) == 0) {
    return(structure(
      list(sat = matrix(0, 1, 1), empty = TRUE),
      class = "NumMatrix"
    ))
  }
  rows <- nrow(img); cols <- ncol(img)
  sat <- matrix(0, nrow = rows + 1, ncol = cols + 1)
  for (i in seq_len(rows)) {
    for (j in seq_len(cols)) {
      sat[i + 1, j + 1] <-
        sat[i, j + 1] + sat[i + 1, j] - sat[i, j] + img[i, j]
    }
  }
  structure(list(sat = sat, empty = FALSE), class = "NumMatrix")
}

#' Sum of a rectangular region in the NumMatrix (0-based inclusive indexing).
sumRegion <- function(obj, row1, col1, row2, col2) {
  if (!inherits(obj, "NumMatrix")) stop("obj must be a NumMatrix")
  if (obj$empty) return(0)
  sat <- obj$sat
  max.r <- dim(sat)[1]; max.c <- dim(sat)[2]
  sat[min(row2 + 2, max.r), min(col2 + 2, max.c)] -
    sat[min(row2 + 2, max.r), col1 + 1] -
    sat[row1 + 1,             min(col2 + 2, max.c)] +
    sat[row1 + 1,             col1 + 1]
}
