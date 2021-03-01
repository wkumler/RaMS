
test_that("error if unable to find files", {
  expect_error(grabMSdata(files = rep(tempfile(), 7), verbosity="very"))
})

test_that("error if no files", {
  expect_error(grabMSdata(files = character(), verbosity="very"))
})

test_that("error if weird files", {
  file.create("blah.txt")
  expect_error(grabMSdata(files = "blah.txt", verbosity="very"))
  file.remove("blah.txt")
})



test_that("checkOutputQuality detects things", {
  grab_what <- "everything"
  output_data <- list(MS1=data.table(runif(100)))
  expect_error(checkOutputQuality(output_data, grab_what))

  grab_what <- c("MS1", "MS2")
  output_data <- list(MS1=data.table(runif(100)), MS2=data.table())
  expect_message(checkOutputQuality(output_data, grab_what))
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