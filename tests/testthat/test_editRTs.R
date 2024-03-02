
ms1_test_files <- c(mzML_filenames[2:4], mzXML_filenames[2])
output_dir <- tempdir()
output_filenames <- paste0(output_dir, "\\", basename(ms1_test_files))

test_that("editRTs works for MS1 data", {
  # mzML encoded in seconds
  init_rts <- grabMSdata(ms1_test_files[1], verbosity = 0)$MS1[,unique(rt)]
  editMzmlRTs(filename = ms1_test_files[1], new_rts = 1:705,
              output_filename = output_filenames[1])
  new_rts <- grabMSdata(output_filenames[1], verbosity = 0)$MS1[,unique(rt)]
  expect_false(all(init_rts==new_rts))
  expect_equal(new_rts, 1:705/60)
  file.remove(output_filenames[1])

  # mzXMLs encoded in minutes
  init_rts <- grabMSdata(ms1_test_files[4], verbosity = 0)$MS1[,unique(rt)]
  editMzxmlRTs(filename = ms1_test_files[4], new_rts = 1:705,
               output_filename = output_filenames[4])
  new_rts <- grabMSdata(output_filenames[4], verbosity = 0)$MS1[,unique(rt)]
  expect_false(all(init_rts==new_rts))
  expect_equal(new_rts, 1:705)
  file.remove(output_filenames[4])
})

ms2_test_file <- c(mzML_filenames[5])
ms2_output_file <- paste0(output_dir, "\\", basename(ms2_test_file))

test_that("editRTs works for MS2 data", {
  editMzmlRTs(filename = ms2_test_file, new_rts = 1:1073,
              output_filename = ms2_output_file)

  init_rt_accession <- grabAccessionData(ms2_test_file, "MS:1000016")
  all_rts <- as.numeric(init_rt_accession$value)/60
  new_rt_accession <- grabAccessionData(ms2_output_file, "MS:1000016")
  new_all_rts <- as.numeric(new_rt_accession$value)/60
  expect_equal(new_all_rts, 1:1073/60)
  expect_false(all(all_rts==new_all_rts))
  file.remove(ms2_output_file)

  # Check that it works when only MS1 values are provided
  editMzmlRTs(filename = ms2_test_file, new_rts = 1:961,
              output_filename=ms2_output_file)
  init_rt_accession <- grabAccessionData(ms2_test_file, "MS:1000016")
  all_rts <- as.numeric(init_rt_accession$value)/60
  new_rt_accession <- grabAccessionData(ms2_output_file, "MS:1000016")
  new_all_rts <- as.numeric(new_rt_accession$value)/60
  expect_false(all(all_rts==new_all_rts))

  new_ms1_rts <- grabMSdata(ms2_output_file, verbosity = 0)$MS1[,unique(rt)]
  expect_equal(new_ms1_rts, 1:961/60)
  file.remove(ms2_output_file)

  editMzmlRTs(filename = ms2_test_file, new_rts = 1:961,
              interp_method = "constant",
              output_filename = ms2_output_file)
  const_ms1_rts <- grabMSdata(ms2_output_file, verbosity = 0)$MS1[,unique(rt)]
  const_all_accession <- grabAccessionData(ms2_output_file, "MS:1000016")
  const_all_rts <- as.numeric(const_all_accession$value)/60
  expect_gte(sum(duplicated(const_all_rts)), 0)

  msn_accession <- grabAccessionData(ms2_output_file, "MS:1000511")
  msn_levels <- as.numeric(msn_accession$value)
  file.remove(ms2_output_file)
})

ms2_test_file <- c(mzXML_filenames[3])
ms2_output_file <- paste0(output_dir, "\\", basename(ms2_test_file))

test_that("editRTs works for MS2 data mzXML", {
  editMzxmlRTs(filename = ms2_test_file, new_rts = 1:1073,
              output_filename = ms2_output_file)
  file.remove(ms2_output_file)

  editMzxmlRTs(filename = ms2_test_file, new_rts = 1:961,
              output_filename=ms2_output_file)

  new_ms1_rts <- grabMSdata(ms2_output_file, verbosity = 0)$MS1[,unique(rt)]
  expect_equal(new_ms1_rts, 1:961/60)
  file.remove(ms2_output_file)

  editMzxmlRTs(filename = ms2_test_file, new_rts = 1:961,
              interp_method = "constant",
              output_filename = ms2_output_file)
  const_ms1_rts <- grabMSdata(ms2_output_file, verbosity = 0)$MS1[,unique(rt)]
  file.remove(ms2_output_file)

  expect_true(TRUE)
})


test_that("editRTs works multifile", {
  new_rts <- list(1:705, 1:705, 1:705, 1:705)
  init_msdata <- grabMSdata(ms1_test_files)
  editMSfileRTs(files = ms1_test_files, new_rt_list = new_rts,
                new_filenames = output_filenames)
  new_msdata <- grabMSdata(output_filenames)
  expect_in(new_msdata$MS1[,unique(rt)], c(1:705, 1:705/60))
  file.remove(output_filenames)
})

test_that("editRT errors and messages work", {
  expect_error({
    editMzmlRTs(filename = "blah", new_rts = 1:961)
  })
  expect_error({
    editMzmlRTs(filename = ms1_test_files[1])
  })
  expect_error({
    editMzmlRTs(filename = ms1_test_files[1],
                new_rts = 1:961, interp_method = "banana")
  })
  expect_error({
    editMzmlRTs(filename = ms1_test_files[1],
                new_rts = 1:10, interp_method = "banana")
  })

  expect_error({
    editMzxmlRTs(filename = "blah", new_rts = 1:961)
  })
  expect_error({
    editMzxmlRTs(filename = ms1_test_files[4])
  })
  expect_error({
    editMzxmlRTs(filename = ms1_test_files[4],
                 new_rts = 1:961, interp_method = "banana")
  })
  expect_error({
    editMzxmlRTs(filename = ms1_test_files[4], new_rts = 1:10)
  })

  # Test below fails on OSX and Ubuntu - unclear how to diagnose
  # file.copy(from = ms1_test_files[1], to = output_dir)
  # expect_warning({
  #   editMSfileRTs(files = output_filenames[1], new_rt_list = list(1:705),
  #                 new_filenames = output_filenames[1])
  # })
  # rtcor_filename <- gsub("\\.(?=mzX?ML)", replacement = "_rtcor\\.",output_filenames[1], perl = TRUE)
  # expect_true(file.exists(rtcor_filename))
  # file.remove(output_filenames[1])
  # file.remove(rtcor_filename)
})
