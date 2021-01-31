
mz_bounds <- pmppm(EIC_mz, ppm = EIC_ppm)

test_that("EIC smaller than all_data", {
  expect_lte(nrow(mzML_EIC_data$EIC), nrow(mzML_data$MS1))
  expect_lte(nrow(mzML_EIC_data$EIC_MS2), nrow(mzML_data$MS2))
})

test_that("EIC extracted expected bounds", {
  expect_gte(min(mzML_EIC_data$EIC$rt),  min(mzML_data$MS1$rt))
  expect_lte(max(mzML_EIC_data$EIC$rt),  max(mzML_data$MS1$rt))

  expect_true(all(mzML_EIC_data$EIC$mz > min(mz_bounds)))
  expect_true(all(mzML_EIC_data$EIC$mz < max(mz_bounds)))
})

# For mzXML
test_that("EIC mzXML smaller than all_data", {
  expect_lte(nrow(mzXML_EIC_data$EIC), nrow(mzXML_data$MS1))
  expect_lte(nrow(mzXML_EIC_data$EIC_MS2), nrow(mzXML_data$MS2))
})

test_that("EIC mzXML extracted expected bounds", {
  expect_gte(min(mzXML_EIC_data$EIC$rt),  min(mzXML_data$MS1$rt))
  expect_lte(max(mzXML_EIC_data$EIC$rt),  max(mzXML_data$MS1$rt))

  expect_true(all(mzXML_EIC_data$EIC$mz > min(mz_bounds)))
  expect_true(all(mzXML_EIC_data$EIC$mz < max(mz_bounds)))
})

