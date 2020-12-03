test_that("mzML_MS1 matches mzXML_MS1", {
  expect_equal(mzML_data$MS1$rt, mzXML_data$MS1$rt, tolerance = 0.001)
  expect_equal(mzML_data$MS1$mz, mzXML_data$MS1$mz)
  expect_equal(mzML_data$MS1$int, mzXML_data$MS1$int, tolerance = 1)

  mzML_filenames <- gsub("mzML", "", unique(mzML_data$MS1$filename))
  mzXML_filenames <- gsub("mzXML", "", unique(mzXML_data$MS1$filename))
  expect_equal(mzML_filenames, mzXML_filenames)
})

test_that("mzML_MS2 matches mzXML_MS2", {
  expect_equal(mzML_data$MS2$rt, mzXML_data$MS2$rt, tolerance = 0.001)
  expect_equal(mzML_data$MS2$mz, mzXML_data$MS2$mz)
  expect_equal(mzML_data$MS2$int, mzXML_data$MS2$int, tolerance = 0.1)

  mzML_filenames <- gsub("mzML", "", unique(mzML_data$MS2$filename))
  mzXML_filenames <- gsub("mzXML", "", unique(mzXML_data$MS2$filename))
  expect_equal(mzML_filenames, mzXML_filenames)
})

test_that("mzML_BPC matches mzXML_BPC", {
  expect_equal(mzML_data$BPC$rt, mzXML_data$BPC$rt, tolerance = 0.001)
  expect_equal(mzML_data$BPC$int, mzXML_data$BPC$int, tolerance = 0.1)

  clean_mzML_filenames <- gsub("mzML", "", unique(mzML_data$BPC$filename))
  clean_mzXML_filenames <- gsub("mzXML", "", unique(mzXML_data$BPC$filename))
  expect_equal(clean_mzML_filenames, clean_mzXML_filenames)
})

test_that("mzML_TIC matches mzXML_TIC", {
  expect_equal(mzML_data$TIC$rt, mzXML_data$TIC$rt, tolerance = 0.001)
  expect_equal(mzML_data$TIC$int, mzXML_data$TIC$int, tolerance = 0.1)

  mzML_filenames <- gsub("mzML", "", unique(mzML_data$TIC$filename))
  mzXML_filenames <- gsub("mzXML", "", unique(mzXML_data$TIC$filename))
  expect_equal(mzML_filenames, mzXML_filenames)
})

test_that("mzML_EIC matches mzXML_EIC", {
  expect_equal(mzML_data$EIC$rt, mzXML_data$EIC$rt, tolerance = 0.001)
  expect_equal(mzML_data$EIC$int, mzXML_data$EIC$int, tolerance = 0.1)

  mzML_filenames <- gsub("mzML", "", unique(mzML_data$EIC$filename))
  mzXML_filenames <- gsub("mzXML", "", unique(mzXML_data$EIC$filename))
  expect_equal(mzML_filenames, mzXML_filenames)
})

test_that("mzML_EIC_MS2 matches mzXML_EIC_MS2", {
  expect_equal(mzML_data$EIC_MS2$rt, mzXML_data$EIC_MS2$rt, tolerance = 0.001)
  expect_equal(mzML_data$EIC_MS2$int, mzXML_data$EIC_MS2$int, tolerance = 0.1)

  mzML_filenames <- gsub("mzML", "", unique(mzML_data$EIC_MS2$filename))
  mzXML_filenames <- gsub("mzXML", "", unique(mzXML_data$EIC_MS2$filename))
  expect_equal(mzML_filenames, mzXML_filenames)
})
