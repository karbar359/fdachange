# Golden values captured from the original (audited) thesis implementation
# by data-raw/capture_golden_values.R at full double precision.
.golden <- local({
  df <- read.csv(test_path("testdata", "golden_values.csv"),
                 stringsAsFactors = FALSE)
  stats::setNames(as.numeric(df$value), df$id)
})

golden <- function(id) {
  v <- .golden[[id]]
  if (is.null(v) || is.na(v)) stop("no golden value for id: ", id)
  v
}

expect_golden <- function(id, value, tolerance = 1e-12) {
  expect_equal(unname(value), golden(id), tolerance = tolerance, label = id)
}
