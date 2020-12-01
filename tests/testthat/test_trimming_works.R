test_rtrange <- c(5, 6)
mzML_trimmed <- grabMSdata(files = mzML_filenames, grab_what = c("MS1", "MS2", "TIC", "BPC"),
                           verbosity = "none", rtrange = test_rtrange)

test_that("trim edges within data", {
  expect_gte(min(test_rtrange), min(mzML_data$MS1$rt))
  expect_lte(max(test_rtrange), max(mzML_data$MS1$rt))
})

test_that("trimming reduces size", {
  expect_gte(nrow(mzML_data$MS1), nrow(mzML_trimmed$MS1))
  expect_gte(nrow(mzML_data$MS2), nrow(mzML_trimmed$MS2))
})

test_that("edges properly trimmed", {
  expect_lte(min(test_rtrange), min(mzML_trimmed$MS1$rt))
  expect_gte(max(test_rtrange), max(mzML_trimmed$MS1$rt))
})

test_rtrange <- c(5, 5)
mzML_trimmed <- grabMSdata(files = mzML_filenames, grab_what = c("MS1", "MS2", "TIC", "BPC"),
                           verbosity = "none", rtrange = test_rtrange)
test_that("zero-length trim returns nothing", {
  expect_equal(mean(sapply(mzML_trimmed, nrow)), 0)
})
