test_that("telco_far1 has 11 clusters of 50 x 28 non-negative curves", {
  expect_named(telco_far1, paste0("cl", 1:11))
  for (Y in telco_far1) {
    expect_true(is.matrix(Y) && is.double(Y))
    expect_equal(dim(Y), c(50L, 28L))
    expect_false(anyNA(Y))
    expect_true(all(Y >= 0))
    expect_named(attributes(Y), c("dim", "dimnames", "change_months",
                                  "multiplier", "jump", "outliers",
                                  "outlier_months"))
    expect_length(attr(Y, "multiplier"), 50L)
    expect_length(attr(Y, "jump"), 50L)
  }
})

test_that("telco_far1 designs match the documentation", {
  expect_equal(attr(telco_far1$cl2, "change_months"), 37:46)
  expect_equal(attr(telco_far1$cl3, "change_months"), 37:49)
  expect_equal(attr(telco_far1$cl4, "change_months"), 37:50)
  for (nm in c("cl1", paste0("cl", 5:11))) {
    expect_length(attr(telco_far1[[nm]], "change_months"), 0L)
    expect_true(all(attr(telco_far1[[nm]], "jump") == 0))
  }
  expect_equal(attr(telco_far1$cl3, "outlier_months"), c(37L, 45L, 49L))
  expect_equal(attr(telco_far1$cl9, "outlier_months"), c(12L, 16L, 30L, 32L))
  expect_equal(nrow(attr(telco_far1$cl9, "outliers")), 64L)
  # the jump raises the monthly mean inside cl2's window
  m <- rowMeans(telco_far1$cl2)
  expect_gt(mean(m[37:46]), mean(m[-(37:46)]))
})
