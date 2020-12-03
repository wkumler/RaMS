
test_that("trim edges within data", {
  expect_gte(min(test_rtrange), min(mzML_data$MS1$rt))
  expect_lte(max(test_rtrange), max(mzML_data$MS1$rt))
})

test_that("trimming reduces size", {
  expect_gte(nrow(mzML_data$MS1), nrow(mzML_trimmed$MS1))
})

test_that("edges properly trimmed", {
  expect_lte(min(test_rtrange), min(mzML_trimmed$MS1$rt))
  expect_gte(max(test_rtrange), max(mzML_trimmed$MS1$rt))
})

test_that("zero-length trim returns nothing", {
  expect_true(all(sapply(mzML_zeroed, nrow)==0))
})

