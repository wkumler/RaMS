test_that("everything grabs everything", {
  expect_identical(mzML_everything, mzML_data)
})

test_that("everything got grabbed", {
  expect_setequal(names(mzML_data), c("MS1", "MS2", "metadata", "TIC", "BPC"))
})

test_that("mzML MS1 data looks normal", {
  # Check that there's data there at all
  expect_gt(nrow(mzML_data$MS1), 0)

  # Check that only expected columns exist
  expect_setequal(names(mzML_data$MS1), c("rt", "mz", "int", "filename"))

  # Check on column types
  expect_type(mzML_data$MS1$rt, "double")
  expect_type(mzML_data$MS1$mz, "double")
  expect_type(mzML_data$MS1$int, "double")
  expect_type(mzML_data$MS1$filename, "character")

  # Check that the expected files were loaded
  expect_length(unique(mzML_data$MS1$filename), 2)
  expect_identical(unique(mzML_data$MS1$filename), basename(mzML_filenames[1:2]))
})

test_that("mzML MS2 data looks normal", {
  # Check that there's data there at all
  expect_gt(nrow(mzML_data$MS2), 0)

  # Check that only expected columns exist
  expect_setequal(names(mzML_data$MS2), c("rt", "fragmz", "premz", "voltage",
                                          "int", "filename"))

  # Check on column types
  expect_type(mzML_data$MS2$rt, "double")
  expect_type(mzML_data$MS2$fragmz, "double")
  expect_type(mzML_data$MS2$premz, "double")
  expect_type(mzML_data$MS2$int, "double")
  expect_type(mzML_data$MS2$voltage, "integer")
  expect_identical(unique(mzML_data$MS2$voltage), 35L)
  expect_type(mzML_data$MS2$filename, "character")

  # Check that the expected files were loaded
  expect_length(unique(mzML_data$MS2$filename), 1)
  expect_identical(unique(mzML_data$MS2$filename), basename(mzML_filenames[1]))
})

test_that("mzML BPC/TIC looks normal", {
  # Check that there's data there at all
  expect_gt(nrow(mzML_data$BPC), 0)

  # Check that only expected columns exist
  expect_setequal(names(mzML_data$BPC), c("rt", "int", "filename"))

  # Check column types
  expect_type(mzML_data$BPC$rt, "double")
  expect_type(mzML_data$BPC$int, "double")
  expect_type(mzML_data$BPC$filename, "character")



  # Check that there's data there at all
  expect_gt(nrow(mzML_data$TIC), 0)

  # Check that only expected columns exist
  expect_setequal(names(mzML_data$TIC), c("rt", "int", "filename"))

  # Check column types
  expect_type(mzML_data$TIC$rt, "double")
  expect_type(mzML_data$TIC$int, "double")
  expect_type(mzML_data$TIC$filename, "character")



  # Check that TIC <= BPC
  expect_true(all(mzML_data$BPC$int<=mzML_data$TIC$int))
})

test_that("multifile works", {
  expect_length(unique(mzML_multi_data$MS1$filename), 3)
  expect_length(unique(mzML_multi_data$TIC$filename), 3)
})
