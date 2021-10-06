
test_that("prefilter didn't ruin things", {
  expect_gt(nrow(mzML_data$MS1), 0)
  expect_gt(nrow(mzML_data$MS2), 0)
  expect_gt(nrow(mzML_data$BPC), 0)
  expect_gt(nrow(mzML_data$TIC), 0)

  expect_gt(nrow(mzXML_data$MS1), 0)
  expect_gt(nrow(mzXML_data$MS2), 0)
  expect_gt(nrow(mzXML_data$BPC), 0)
  expect_gt(nrow(mzXML_data$TIC), 0)
})

test_that("negative prefilter doesn't change anything", {
  expect_equal(nrow(mzML_data$MS1), nrow(mzML_nofilter$MS1))
})

test_that("prefilter removes data points", {
  expect_gte(nrow(mzML_data$MS1), nrow(mzML_majorfilter$MS1))
  expect_gte(min(mzML_majorfilter$MS1$int), major_prefilter)
})

test_that("prefilter removes expected data", {
  expect_equal(mzML_data$MS1[int>major_prefilter], mzML_majorfilter$MS1)
})

