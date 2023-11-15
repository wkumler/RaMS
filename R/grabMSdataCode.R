# Welcome to RaMS!

# grabMSdata ----

#' Grab mass-spectrometry data from file(s)
#'
#' The main `RaMS` function. This function accepts a list of the files that will
#' be read into R's working memory and returns a list of `data.table`s
#' containing the requested information. What information is requested is
#' determined by the `grab_what` argument, which can include MS1, MS2, BPC, TIC,
#' or metadata information. This function serves as a wrapper around both
#' `grabMzmlData` and `grabMzxmlData` and handles multiple files, but those two
#' have also been exposed to the user in case super-simple handling is desired.
#' Retention times are reported in minutes, and will be converted automatically
#' if they are encoded in seconds.
#'
#' @param files A character vector of filenames to read into R's memory. Both
#'   absolute and relative paths are acceptable.
#' @param grab_what What data should be read from the file? Options include
#'   "MS1" for data only from the first spectrometer, "MS2" for fragmentation
#'   data, "BPC" for rapid access to the base peak chromatogram, "TIC" for rapid
#'   access to the total ion chromatogram, "DAD" for DAD (UV) data, and "chroms"
#'   for precompiled chromatogram data (especially useful for MRM but often
#'   contains BPC/TIC in other files). Metadata can be accessed with "metadata",
#'   which provides information about the instrument and time the file was run.
#'   These options can be combined (i.e. `grab_data=c("MS1", "MS2", "BPC")`) or
#'   this argument can be set to "everything" to extract all of the above.
#'   Options "EIC" and "EIC_MS2" are useful when working with files whose total
#'   size exceeds working memory - it first extracts all relevant MS1 and MS2
#'   data, respectively, then discards data outside of the mass range(s)
#'   calculated from the provided mz and ppm. The default, "everything",
#'   includes all MS1, MS2, BPC, TIC, and metadata.
#' @param verbosity Three levels of processing output to the R console are
#'   available, with increasing verbosity corresponding to higher integers. A
#'   verbosity of zero means that no output will be produced, useful when
#'   wrapping within larger functions. A verbosity of 1 will produce a progress
#'   bar using base R's txtProgressBar function. A verbosity of 2 or higher will
#'   produce timing output for each individual file read in. The default, NULL,
#'   will select between 1 and 2 depending on the number of files being read: if
#'   a single file, verbosity is set to 2; if multiple files, verbosity is set
#'   to 1.
#' @param mz A vector of the mass-to-charge ratio for compounds of interest.
#'   Only used when combined with `grab_what = "EIC"` (see above). Multiple
#'   masses can be provided.
#' @param ppm A single number corresponding to the mass accuracy (in parts per
#'   million) of the instrument on which the data was collected. Only used when
#'   combined with `grab_what = "EIC"` (see above).
#' @param rtrange Only available when parsing mzML files. A vector of length 2
#'   containing an upper and lower bound on retention times of interest.
#'   Providing a range here can speed up load times (although not enormously, as
#'   the entire file must still be read) and reduce the final object's size.
#' @param prefilter A single number corresponding to the minimum intensity of
#'   interest in the MS1 data. Data points with intensities below this threshold
#'   will be silently dropped, which can dramatically reduce the size of the
#'   final object. Currently only works with MS1 data, but could be expanded
#'   easily to handle more.
#'
#' @return A list of `data.table`s, each named after the arguments requested in
#'   grab_what. E.g. $MS1 contains MS1 information, $MS2 contains fragmentation
#'   info, etc. MS1 data has four columns: retention time (rt), mass-to-charge
#'   (mz), intensity (int), and filename. MS2 data has six: retention time (rt),
#'   precursor m/z (premz), fragment m/z (fragmz), fragment intensity (int),
#'   collision energy (voltage), and filename. Data requested that does not
#'   exist in the provided files (such as MS2 data requested from MS1-only
#'   files) will return an empty (length zero) data.table. The data.tables
#'   extracted from each of the individual files are collected into one large
#'   table using data.table's `rbindlist`. $metadata is a little weirder because
#'   the metadata doesn't fit neatly into a tidy format but things are hopefully
#'   named helpfully. $chroms was added in v1.3 and contains 7 columns:
#'   chromatogram type (usually TIC, BPC or SRM info), chromatogram index,
#'   target mz, product mz, retention time (rt), and intensity (int). $DAD was
#'   also added in v1.3 and contains has three columns: retention time (rt),
#'   wavelength (lambda),and intensity (int). Data requested that does not exist
#'   in the provided files (such as MS2 data requested from MS1-only files) will
#'   return an empty (zero-row) data.table.
#'
#' @export
#'
#' @examples
#' library(RaMS)
#' # Extract MS1 data from a couple files
#' sample_dir <- system.file("extdata", package = "RaMS")
#' sample_files <- list.files(sample_dir, full.names=TRUE)
#' multifile_data <- grabMSdata(sample_files[c(3, 5, 6)], grab_what="MS1")
#'
#' # "Stream" data from the internet (i.e. Metabolights)
#' \dontrun{
#' access_url <- "https://www.ebi.ac.uk/metabolights/MTBLS703/files"
#'
#' # URL below obtained by right-clicking site download button and copying
#' # link address
#' sample_url <- paste0("https://www.ebi.ac.uk/metabolights/ws/studies/",
#'                      "MTBLS703/download/acefcd61-a634-4f35-9c3c-c572",
#'                      "ade5acf3?file=161024_Smp_LB12HL_AB_pos.mzXML")
#' file_data <- grabMSdata(sample_url, grab_what="everything", verbosity=2)
#' }
grabMSdata <- function(files, grab_what="everything", verbosity=NULL,
                       mz=NULL, ppm=NULL, rtrange=NULL, prefilter=-1){
  # Check that files were provided
  if(!length(files)>0)stop("No files provided")

  # Check that grab_what is one of the approved options
  good_grabs <- c("MS1", "MS2", "EIC", "EIC_MS2", "everything", "metadata",
                  "BPC", "TIC", "chroms", "DAD")
  if(any(!grab_what%in%good_grabs)){
    bad_grabs <- paste(grab_what[!grab_what%in%good_grabs], collapse = ", ")
    stop(paste0("`grab_what = ", bad_grabs, "` is not currently supported"))
  }

  # Handle tmzMLs first
  tmzml_check <- grepl("\\.tmzML", files)
  if(any(tmzml_check)){
    msdata_con_object <- msdata_connection_constructor(
      tmzml_check, files, grab_what, mz, ppm, rtrange, prefilter, verbosity
    )
    return(msdata_con_object)
  }

  # Handle mzMLs and mzXMLs below
  # Check for non-mzXML or non-mzML files
  checkQuickMLs(files)

  # Add sanity check for EIC extraction
  if(!is.null(mz) & !any(c("EIC", "EIC_MS2")%in%grab_what)){
    warning(paste0('Argument "mz" should be used with grab_what = "EIC" or',
                   '"EIC_MS2" and will be ignored in the current call'))
  }
  if(!is.null(ppm) & !any(c("EIC", "EIC_MS2")%in%grab_what)){
    warning(paste0('Argument "mz" should be used with grab_what = "EIC" or',
                   '"EIC_MS2" and will be ignored in the current call'))
  }

  # Handle null verbosity flag with intelligent defaults
  if(is.null(verbosity)){
    verbosity <- ifelse(length(files)==1, 2, 1)
  }

  # Define outer control loop so multiple files can be read in simultaneously
  all_file_data <- list()
  if(verbosity>0){
    if(length(files)>=2){
      pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    }
    start_time <- Sys.time()
  }

  for(i in seq_along(files)){
    filename <- files[i]

    if(grepl("\\.mzML", basename(filename), ignore.case = TRUE)){
      out_data <- grabMzmlData(filename = filename, grab_what = grab_what,
                               verbosity = verbosity, mz = mz, ppm = ppm,
                               rtrange = rtrange, prefilter = prefilter)
    } else if(grepl("\\.mzXML", basename(filename), ignore.case = TRUE)){
      out_data <- grabMzxmlData(filename = filename, grab_what = grab_what,
                                verbosity = verbosity, mz = mz, ppm = ppm,
                                rtrange = rtrange, prefilter = prefilter)
    } else {
      stop(paste("Unable to determine file type for", filename))
    }
    out_data_filenamed <- mapply(function(dt_i, fname_i){
      if(nrow(dt_i)){
        dt_i$filename <- fname_i
      } else {
        dt_i$filename <- character()
      }
      dt_i
    }, dt_i=out_data, fname_i=basename(filename), SIMPLIFY = FALSE)
    all_file_data[[i]] <- out_data_filenamed
    names(all_file_data)[[i]] <- basename(filename)
    if(verbosity>0 & length(files)>=2){
      setTxtProgressBar(pb, i)
    }
  }
  if(verbosity>0){
    if(length(files)>=2){
      close(pb)
    }
    if(verbosity>1){
      cat("Binding files together into single data.table\n")
    }
  }

  # Bind all the similar pieces together (e.g. stack MS1 from different files)
  # all_file_data_output <- Reduce(function(x,y) Map(rbind, x, y, fill=TRUE), all_file_data)
  all_file_data_output <- lapply(seq_along(all_file_data[[1]]), function(sub_index){
    rbindlist(lapply(all_file_data, `[[`, sub_index))
  })
  names(all_file_data_output) <- names(all_file_data[[1]])

  invisible(checkOutputQuality(
    output_data = all_file_data_output, grab_what = grab_what
  ))

  if(verbosity>0){
    time_total <- round(difftime(Sys.time(), start_time), digits = 2)
    cat("Total time:", time_total, units(time_total), "\n")
  }

  all_file_data_output
}

# checkOutputQuality ----

#' Check that the output data is properly formatted.
#'
#' This function checks that data produced by repeated calls to the
#' `grabMzmlData()` and `grabMzxmlData()` functions is formatted properly before
#' it's provided to the user. It checks that all of the requested data has been
#' obtained and warns if data is found to be empty, misnamed, or has columns of
#' the wrong type.
#'
#' @param output_data The collected data resulting from repeated calls to
#'   `grabMzmlData()`, after being bound together.
#' @param grab_what The names of the data requested by the user.
#'
#' @return NULL (invisibly). The goal of this function is its side effects, i.e.
#'   throwing errors and providing info when the files are not found.
#'
checkOutputQuality <- function(output_data, grab_what){
  if("everything"%in%grab_what){
    grab_what <- c("MS1", "MS2", "BPC", "TIC", "metadata")
  }
  missing_data <- !grab_what%in%names(output_data)
  if(any(missing_data)){
    warning(paste("Not all data collected; missing",
                  paste(grab_what[missing_data], collapse = ", ")))
  }

  missing_data <- !as.logical(sapply(output_data, length))
  if(any(missing_data)){
    message(paste("Some bits seem to be empty:",
                  paste(names(output_data)[missing_data], collapse = ", ")))
  }

  zero_rows <- !as.logical(sapply(output_data, nrow))
  if(any(missing_data)){
    message(paste("Some data seem to be empty:",
                  paste(names(output_data)[missing_data], collapse = ", ")))
  }

  misnamed <- mapply(function(nms, dt){
    if(nms=="BPC"|nms=="TIC"){
      proper_names <- c("rt", "int", "filename")
    } else if(nms=="MS1"){
      proper_names <- c("rt", "mz", "int", "filename")
    } else if(nms=="MS2"){
      proper_names <- c("rt", "premz", "fragmz", "int", "voltage", "filename")
    } else if (nms=="DAD"){
      proper_names <- c("rt", "lambda", "int", "filename")
    } else if (nms=="EIC"){
      proper_names <- c("rt", "mz", "int", "filename")
    } else if (nms=="EIC_MS2"){
      proper_names <- c("rt", "premz", "fragmz", "int", "voltage", "filename")
    } else if(nms=="chroms"){
      proper_names <- c("chrom_type", "chrom_index", "target_mz", "product_mz",
                        "rt", "int")
    } else if(nms=="metadata"){
      return(FALSE) # Because we don't know what names metadata may have
    } else {
      message(paste0("Unexpected data subset name in ",
                     "checkQualityOutput, may be malformed"))
      proper_names <- nms
    }
    !all(proper_names%in%names(dt))&all(names(dt)%in%proper_names)
  }, nms=names(output_data), dt=output_data)
  if(any(misnamed)){
    message(paste("Some data seem to have incorrect names:",
                  paste(names(output_data)[misnamed], collapse = ", ")))
  }

  wrong_coltype <- sapply(output_data, function(dt){
    col_classes <- sapply(dt, class)
    proper_class <- c(rt="numeric", mz="numeric", int="numeric",
                      filename="character", premz="numeric", fragmz="numeric",
                      voltage="integer")
    which(proper_class[names(col_classes)] != col_classes)
  })
  if(any(sapply(wrong_coltype, length))){
    message("Some data seem to have incorrect column types:")
    invisible(mapply(function(incorrect_coltypes, incorrect_names){
      if(length(incorrect_coltypes)){
        message(paste(incorrect_names, names(incorrect_coltypes), "\n"))
      }
    }, wrong_coltype, names(wrong_coltype)))
  }
  return(invisible(NULL))
}

# grabAccessionData ----
#' Get arbitrary metadata from an mzML file by accession number
#'
#' @param filename The name of the file for which metadata is requested. Both
#'   absolute and relative paths are acceptable.
#' @param accession_number The HUPO-PSI accession number for the metadata to
#' be extracted. These accession numbers are typically of the form MS:#######
#' and the full list can be found and searched at
#' https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo.
#'
#' @return A data frame with the name and value of the parameter requested,
#' as deduced from the XML tag attributes corresponding to the accession
#' number.
#' @export
#'
#' @examples
#' library(RaMS)
#' sample_dir <- system.file("extdata", package = "RaMS")
#' sample_file <- list.files(sample_dir, full.names=TRUE)[3]
#' # Get ion injection time
#' iit_df <- grabAccessionData(sample_file, "MS:1000927")
#' # Manually create TIC
#' int_df <- grabAccessionData(sample_file, "MS:1000285")
#' rt_df <- grabAccessionData(sample_file, "MS:1000016")
#' tic <- data.frame(rt=rt_df$value, int=int_df$value)
#' plot(tic$rt, tic$int, type = "l")
grabAccessionData <- function(filename, accession_number){
  xml_data <- xml2::read_xml(filename)

  checkNamespace(xml_data)
  arb_xpath <- paste0('//d1:cvParam[@accession="', accession_number, '"]')
  arb_nodes <- xml2::xml_find_all(xml_data, arb_xpath)
  out_df <- data.frame(name=xml2::xml_attr(arb_nodes, "name"),
                       value=xml2::xml_attr(arb_nodes, "value"))
  if(nrow(out_df)==0){
    warning(paste("Accession number", accession_number, "not found"))
  }
  out_df
}

# msdata_connection_constructor ----
msdata_connection_constructor <- function(tmzml_check, files, grab_what, mz, ppm,
                                          rtrange, prefilter, verbosity){
  if(!all(tmzml_check)){
    stop("At this time, tmzMLs can't be mixed with mzML/mzXMLs")
  }
  if(grab_what=="everything"){
    grab_what <- c("MS1", "MS2")
  }
  if(!all(grab_what%in%c("MS1", "MS2"))){
    stop("At this time, tmzMLs can only be used with MS1 or MS2 data")
  }
  if(!is.null(mz)){
    warning("Argument 'mz' has no function when used with tmzML files, ignoring")
  }
  if(!is.null(ppm)){
    warning("Argument 'ppm' has no function when used with tmzML files, ignoring")
  }
  if(!is.null(rtrange)){
    warning("Argument 'rtrange' has no function when used with tmzML files, ignoring")
  }
  if(prefilter>-1){
    warning("Argument 'prefilter' has no function when used with tmzML files, ignoring")
  }

  # Handle null verbosity flag with intelligent defaults
  if(is.null(verbosity)){
    verbosity <- ifelse(length(files)==1, 0, 1)
  }

  # Check for missing files before creating object
  url_files <- grepl(x = files, "^(http|ftp)")
  file_exists <- file.exists(files[!url_files])
  if(!all(file_exists)){
    stop(paste("Unable to find all files, e.g.\n",
               paste(head(files[!url_files][!file_exists]), collapse = "\n ")))
  }

  # Create a list object to hide connection values and allow
  # RStudio to autocomplete MS1 and MS2
  msdata_con <- vector("list", length = length(grab_what)+1)
  msdata_con[[length(msdata_con)]] <- list(
    files=files,
    grab_what=grab_what,
    verbosity=verbosity
  )
  names(msdata_con) <- c(grab_what, "connection")

  class(msdata_con) <- "msdata_connection"
  return(msdata_con)
}
# Other functions ----

checkNamespace <- function(xml_data){
  if (is.na(xml2::xml_attr(xml_data,"xmlns"))){
    xml2::xml_attr(xml_data,"xmlns") <- "http://psi.hupo.org/ms/mzml"
  }
}

checkFileType <- function(xml_data, node_to_check){
  # Check for mzML node
  # Length works because external pointer has length 2
  if(!length(xml2::xml_find_first(xml_data, paste0("//d1:", node_to_check)))){
    stop(paste0("No ", node_to_check, " node found in this file"))
  }
}

checkProvidedMzPpm <- function(mz, ppm){
  if(is.null(mz)){
    stop("Please provide an m/z value when using grab_what = EIC")
  }
  if(!inherits(mz, "numeric")&&!inherits(mz, "integer")){
    stop("Please provide a numeric m/z value")
  }
  if(any(is.na(mz))){
    stop("It looks like you have an NA among your m/z input")
  }
  if(any(mz<0)){
    stop("m/z must be positive")
  }

  if(is.null(ppm)){
    stop("Please provide a ppm value when using grab_what = EIC")
  }
  if(!inherits(ppm, "numeric")&&!inherits(ppm, "integer")){
    stop("Please provide a numeric ppm value")
  }
  if(ppm<0){
    stop("ppm must be positive")
  }
}

checkRTrange <- function(rtrange){
  if(!is.null(rtrange)){
    if("matrix"%in%class(rtrange)){
      rtrange <- as.vector(rtrange)
    }
    if(length(rtrange)!=2){
      stop("Please provide an rtrange of length 2")
    }
    if(!inherits(rtrange, "numeric")&&!inherits(rtrange, "integer")){
      stop("Please provide a numeric rtrange")
    }
  }
  rtrange
}

checkProvidedPrefilter <- function(prefilter){
  if(!is.numeric(prefilter)){
    warning(paste0("`prefilter` argument should be numeric, but instead it ",
                   "seems to be ", class(prefilter), " with value ",
                   prefilter, ". Replacing with -1..."))

    prefilter <- -1
  }
  if(length(prefilter)>1){
    warning(paste0("`prefilter` argument should be length 1, but instead it ",
                   "seems to be ", length(prefilter), ". Ignoring all extra ",
                   "values..."))
    prefilter <- prefilter[1]
  }
  prefilter
}

checkQuickMLs <- function(files){
  init_warn_val <- getOption("warn")
  options(warn=1)
  ml_check <- grepl("\\.mzx?ml$|\\.mzx?ml\\.gz$", files, ignore.case = TRUE)
  if(!all(ml_check)){
    nonml_files <- files[!ml_check]
    if(length(nonml_files)>5){
      warn_string <- paste0(
        "Not all files appear to be mzML or mzXML: see\n",
        paste(head(nonml_files, 5), collapse = "\n"),
        "\nand ", length(nonml_files)-5, " others"
      )
    } else {
      warn_string <- paste0(
        "Not all files appear to be mzML or mzXML: see\n",
        paste(nonml_files, collapse = "\n")
      )
    }
    warning(warn_string)
  }
  options(warn=init_warn_val)
}

timeReport <- function(last_time, text=NULL){
  time_total <- round(difftime(Sys.time(), last_time), digits = 2)
  cat(time_total, units(time_total), "\n")
  cat(text)
  Sys.time()
}


#' Plus/minus parts per million
#'
#' It shouldn't be hard to translate a point mass into a mass window bounded by
#' spectrometer accuracy.
#'
#' @param mass A length-1 numeric representing the mass of
#' interest for which a mass range is desired.
#' @param ppm The parts-per-million accuracy of the mass spectrometer on which
#' the data was collected.
#'
#' @return A length-2 numeric representing the mass range requested
#'
#' @export
#'
#' @examples
#' pmppm(100, 5)
#' pmppm(1000000, 5)
#' pmppm(118.0865, 2.5)
#' pmppm(892.535313, 10)
pmppm <- function(mass, ppm=4){
  if(length(mass)>1)warning("Vectorized masses are not supported in pmppm")
  if(length(ppm)>1)warning("Vectorized ppms are not supported in pmppm")
  c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))
}

#' Trapezoidal integration of mass-spec retention time / intensity values
#'
#' Performs a trapezoidal Riemann sum to calculate the area under the curve
#' for mass-spectrometry data. Accepts a vector of retention times and the
#' associated intensities and returns the area.
#'
#' @param rts A numeric vector of retention times across an MS peak. Should be
#' monotonically increasing and without duplicates or will throw a warning.
#' @param ints A numeric vector of measured intensities across an MS peak
#' @param baseline A length-1 character vector of either "none" (the default),
#' "square", or "trapezoid".
#'
#' @return A length-1 numeric value representing the area under the curve
#' @export
#'
#' @examples
#' trapz(1:10, 1:10)
#' trapz(1:10, 10:1)
#'
#' trapz(1:10, 11:20)
#' trapz(1:10, 11:20, baseline="square")
#' trapz(1:10, 11:20, baseline="trapezoid")
#'
#' x_vals <- seq(-2, 2, length.out=100)
#' trapz(x_vals, dnorm(x_vals))
#'
#' x_vals <- seq(0, pi/2, length.out=100)
#' trapz(x_vals, cos(x_vals))
trapz <- function(rts, ints, baseline="none"){
  # Checks for argument accuracy
  if(!baseline%in%c("none", "square", "trapezoid")){
    stop("Argument 'baseline' must be one of 'none', 'square', or 'trapezoid'")
  }
  if(length(rts)!=length(ints))stop("rts length must be equal to ints length")
  if(!is.numeric(rts))stop("rts vector must be numeric")
  if(!is.numeric(ints))stop("ints vector must be numeric")
  if(length(rts)<2)return(0)

  # Warn if duplicated RT entries (maybe better as any(duplicated(rts)) ?)
  if(length(rts)!=length(unique(rts)))warning("RTs appear to be duplicated")
  # Warn if RTs do not increase monotonically
  # Test stolen from https://stackoverflow.com/q/13093912
  if(!all(rts == cummax(rts)))warning("RTs do not appear to be ordered")
  # Warn with very large peak widths
  if(length(rts)>1000)warning("Most peak widths are less than 1000 scans wide")
  # Warn with negative intensity values
  if(any(ints<0))warning("Detected intensities less than zero")

  # Calculate the area of each trapezoid (area under curve)
  rt_wedges <- rts[-1]-rts[-length(rts)]
  int_wedges <- (ints[-1]+ints[-length(ints)])/2
  auc <- sum(rt_wedges*int_wedges)

  # Calculate area under the baseline (area under baseline)
  if(baseline=="square"){
    aub <- min(ints)*(max(rts)-min(rts))
  } else if(baseline=="trapezoid"){
    minmax_rts <- c(rts[1], rts[length(rts)])
    minmax_ints <- c(ints[1], ints[length(ints)])
    aub <- trapz(minmax_rts, minmax_ints, baseline="none")
  } else {
    aub <- 0
  }

  return(auc-aub)
}


timeReport <- function(last_time, text=NULL){
  time_total <- round(difftime(Sys.time(), last_time), digits = 2)
  cat(time_total, units(time_total), "\n")
  cat(text)
  Sys.time()
}


# Import area ----

#' @import xml2
#' @import data.table
#' @import utils
#' @importFrom base64enc base64decode base64encode
#' @export
data.table::between
#' @export
data.table::`%between%`
NULL
