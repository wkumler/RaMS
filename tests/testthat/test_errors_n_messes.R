
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
    "MS1", "MS2", "MS3", "MS4"
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

# Minimal documents with one spectrum/scan per polarity given
mini_mzML <- function(pols){
  spectra <- paste0('<spectrum><cvParam name="', pols, '"/></spectrum>', collapse="")
  xml2::read_xml(paste0('<mzML xmlns="http://psi.hupo.org/ms/mzml"><run><spectrumList>',
                        spectra, '</spectrumList></run></mzML>'))
}
mini_mzXML <- function(pols){
  scans <- paste0('<scan polarity="', pols, '"/>', collapse="")
  xml2::read_xml(paste0('<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">',
                        '<msRun>', scans, '</msRun></mzXML>'))
}

test_that("checkMixedPolarity warns only when an mzML has both polarities", {
  expect_silent(checkMixedPolarity(mini_mzML(c("positive scan", "positive scan")), "a.mzML", "mzML"))
  expect_silent(checkMixedPolarity(mini_mzML("negative scan"), "a.mzML", "mzML"))
  expect_silent(checkMixedPolarity(mini_mzML(character(0)), "a.mzML", "mzML"))
  expect_warning(checkMixedPolarity(mini_mzML(c("positive scan", "negative scan")),
                                    "dir/a.mzML", "mzML"),
                 "^a.mzML contains both positive and negative")
})

test_that("checkMixedPolarity warns only when an mzXML has both polarities", {
  expect_silent(checkMixedPolarity(mini_mzXML(c("+", "+")), "a.mzXML", "mzXML"))
  expect_silent(checkMixedPolarity(mini_mzXML("-"), "a.mzXML", "mzXML"))
  expect_silent(checkMixedPolarity(mini_mzXML(character(0)), "a.mzXML", "mzXML"))
  expect_warning(checkMixedPolarity(mini_mzXML(c("+", "-")), "dir/a.mzXML", "mzXML"),
                 "^a.mzXML contains both positive and negative")
})

test_that("checkMixedPolarity uses the XPaths for the given file type", {
  # Each format's polarity markup is invisible to the other format's check
  expect_silent(checkMixedPolarity(mini_mzXML(c("+", "-")), "a.mzXML", "mzML"))
  expect_silent(checkMixedPolarity(mini_mzML(c("positive scan", "negative scan")), "a.mzML", "mzXML"))
})

test_that("checkMixedPolarity flags the mixed-polarity sample files", {
  read_sample <- function(fname){
    xml_data <- xml2::read_xml(system.file("extdata", fname, package="RaMS"))
    if(grepl("mzML", fname))checkNamespace(xml_data)
    xml_data
  }
  expect_warning(checkMixedPolarity(read_sample("S30657.mzML.gz"), "S30657.mzML.gz", "mzML"),
                 "S30657.mzML.gz")
  expect_warning(checkMixedPolarity(read_sample("S30657.mzXML.gz"), "S30657.mzXML.gz", "mzXML"),
                 "S30657.mzXML.gz")
  expect_silent(checkMixedPolarity(read_sample("LB12HL_AB.mzML.gz"), "LB12HL_AB.mzML.gz", "mzML"))
  expect_silent(checkMixedPolarity(read_sample("LB12HL_AB.mzXML.gz"), "LB12HL_AB.mzXML.gz", "mzXML"))
})
