
test_that("error if unable to find files", {
  expect_error(grabMSdata(files = rep(tempfile(fileext = ".mzML"), 7), verbosity=2))
})

test_that("error if no files", {
  expect_error(grabMSdata(files = character(), verbosity=2))
})

test_that("error if grab_what is weird", {
  expect_error(grabMSdata(files = mzML_filenames[2], grab_what = "banana"))
  expect_error(grabMSdata(files = mzML_filenames[2], grab_what = "SM1"))
  expect_error(grabMSdata(files = mzML_filenames[1], grab_what = c(
    "MS1", "MS2", "MS3"
  )))
  expect_error(grabMSdata(files = mzML_filenames[2], grab_what = ""))
})

test_that("checkOutputQuality detects things", {
  grab_what <- "everything"
  output_data <- list(MS1=data.table(runif(100)))
  expect_warning(checkOutputQuality(output_data, grab_what))

  grab_what <- c("MS1", "MS2")
  output_data <- list(MS1=data.table(runif(100)), MS2=data.table())
  expect_message(checkOutputQuality(output_data, grab_what))
})

test_that("mz ppm sanity checks work", {
  expect_warning(grabMSdata(files = mzML_filenames[2],
                            grab_what = "everything",
                            mz=118.0865))
  expect_warning(grabMSdata(files = mzML_filenames[2],
                            grab_what = "everything",
                            ppm=5))
})

test_that("checkProvidedMzPpm detects things", {
  expect_error(checkProvidedMzPpm(mz=NULL))
  expect_error(checkProvidedMzPpm(mz="banana"))
  expect_error(checkProvidedMzPpm(mz=c(100, NA_integer_, 3)))
  expect_error(checkProvidedMzPpm(mz=-3))

  expect_error(checkProvidedMzPpm(mz=100, ppm = NULL))
  expect_error(checkProvidedMzPpm(mz=100, ppm = "banana"))
  expect_error(checkProvidedMzPpm(mz=100, ppm = -3))
})

test_that("checkProvidedPrefilter does things", {
  expect_warning(
    one_out <- checkProvidedPrefilter("banana")
  )
  expect_identical(one_out, -1)

  multi_prefilter <- 1:3
  expect_warning(
    multi_out <- checkProvidedPrefilter(multi_prefilter)
  )
  expect_identical(multi_out, 1L)

  expect_warning(checkProvidedPrefilter(NA))
})

test_that("default verbosity works", {
  expect_output(
    grabMSdata(files = mzML_filenames[2], grab_what = "everything"),
    regexp = "Reading MS1"
  )
})

test_that("checkQuickMLs catches weird file types", {
  expect_warning({
    expect_error(grabMSdata(files = rep(tempfile(), 7)))
  })
  expect_warning({
    expect_error(grabMSdata("desktop.ini"))
  })
})

