# Welcome to RaMS!

# grabMzxmlData ----

#' Get mass-spectrometry data from an mzXML file
#'
#' This function handles the mzXML side of things, reading in files that are
#' written in the mzXML format. Much of the code is similar to the mzXML format,
#' but the xpath handles are different and the mz/int array is encoded
#' simultaneously rather than as two separate entries. This function has been
#' exposed to the user in case per-file optimization (such as peakpicking or
#' additional filtering) is desired before the full data object is returned.
#'
#' @param filename A single filename to read into R's memory. Both absolute and
#'   relative paths are acceptable.
#' @param grab_what What data should be read from the file? Options include
#'   "MS1" for data only from the first spectrometer, "MS2" for fragmentation
#'   data, "BPC" for rapid access to the base peak chromatogram, and "TIC" for
#'   rapid access to the total ion chromatogram. These options can be combined
#'   (i.e. `grab_data=c("MS1", "MS2", "BPC")`) or this argument can be set to
#'   "everything" to extract all of the above. Option "EIC" is useful when
#'   working with files whose total size exceeds working memory - it first
#'   extracts all relevant MS1 and MS2 data, then discards data outside of the
#'   mass range(s) calculated from the provided mz and ppm.
#' @param verbosity Three levels of processing output to the R console are
#'   available, with increasing verbosity corresponding to higher integers. A
#'   verbosity of zero means that no output will be produced, useful when
#'   wrapping within larger functions. A verbosity of 1 will produce a progress
#'   bar using base R's txtProgressBar function. A verbosity of 2 or higher will
#'   produce timing output for each individual file read in.
#' @param mz A vector of the mass-to-charge ratio for compounds of interest.
#'   Only used when combined with `grab_what = "EIC"` (see above). Multiple
#'   masses can be provided.
#' @param ppm A single number corresponding to the mass accuracy (in parts per
#'   million) of the instrument on which the data was collected. Only used when
#'   combined with `grab_what = "EIC"` (see above).
#' @param rtrange Not supported for mzXML data. Only provided here so as to
#'   throw a friendly warning rather than an unexpected error.
#' @param prefilter A single number corresponding to the minimum intensity of
#'   interest in the MS1 data. Data points with intensities below this threshold
#'   will be silently dropped, which can dramatically reduce the size of the
#'   final object. Currently only works with MS1 data, but could be expanded
#'   easily to handle more.
#'
#' @return A list of `data.table`s, each named after the arguments requested in
#'   grab_what. $MS1 contains MS1 information, $MS2 contains fragmentation info,
#'   etc. MS1 data has three columns: retention time (rt), mass-to-charge (mz),
#'   and intensity (int). MS2 data has five: retention time (rt), precursor m/z
#'   (premz), fragment m/z (fragmz), fragment intensity (int), and collision
#'   energy (voltage). Data requested that does not exist in the provided files
#'   (such as MS2 data requested from MS1-only files) will return an empty
#'   (length zero) data.table.
#'
#' @export
#'
#' @examples
#' sample_file <- system.file("extdata", "LB12HL_AB.mzXML.gz", package = "RaMS")
#' file_data <- grabMzxmlData(sample_file, grab_what="MS1")
#' # Extract MS1 data and a base peak chromatogram
#' file_data <- grabMzxmlData(sample_file, grab_what=c("MS1", "BPC"))
#' # Extract EIC for a specific mass
#' file_data <- grabMzxmlData(sample_file, grab_what="EIC", mz=118.0865, ppm=5)
#' # Extract EIC for several masses simultaneously
#' file_data <- grabMzxmlData(sample_file, grab_what="EIC", ppm=5,
#'                            mz=c(118.0865, 146.118104, 189.123918))
#'
#' # Extract MS2 data
#' sample_file <- system.file("extdata", "DDApos_2.mzXML.gz", package = "RaMS")
#' MS2_data <- grabMzxmlData(sample_file, grab_what="MS2")
#'
grabMzxmlData <- function(filename, grab_what, verbosity=0,
                          rtrange=NULL, mz=NULL, ppm=NULL, prefilter=-1){
  if(verbosity>1){
    cat(paste0("\nReading file ", basename(filename), "... "))
    last_time <- Sys.time()
  }
  xml_data <- xml2::read_xml(filename)

  checkFileType(xml_data, "mzXML")
  rtrange <- checkRTrange(rtrange)
  prefilter <- checkProvidedPrefilter(prefilter)

  output_data <- list()

  if("everything"%in%grab_what){
    if(length(setdiff(grab_what, "everything"))&&verbosity>0){
      message(paste("Heads-up: grab_what = `everything` includes",
                    "MS1, MS2, BPC, and TIC data"))
      message("Ignoring additional grabs")
    }
    grab_what <- c("MS1", "MS2", "BPC", "TIC", "metadata")
  }

  if(any(c("MS1", "MS2", "EIC", "EIC_MS2")%in%grab_what)){
    file_metadata <- grabMzxmlEncodingData(xml_data)
  }

  if("MS1"%in%grab_what){
    if(verbosity>1)last_time <- timeReport(last_time, text = "Reading MS1 data...")
    output_data$MS1 <- grabMzxmlMS1(xml_data = xml_data, rtrange=rtrange,
                                    file_metadata = file_metadata, prefilter = prefilter)
  }

  if("MS2"%in%grab_what){
    if(verbosity>1)last_time <- timeReport(last_time, text = "Reading MS2 data...")
    output_data$MS2 <- grabMzxmlMS2(xml_data = xml_data, rtrange=rtrange,
                                    file_metadata = file_metadata)
  }

  if("BPC"%in%grab_what){
    if(verbosity>1)last_time <- timeReport(last_time, text = "Reading BPC...")
    output_data$BPC <- grabMzxmlBPC(xml_data = xml_data, rtrange=rtrange)
  }

  if("TIC"%in%grab_what){
    if(verbosity>1)last_time <- timeReport(last_time, text = "Reading TIC...")
    output_data$TIC <- grabMzxmlBPC(xml_data = xml_data, rtrange=rtrange,
                                    TIC = TRUE)
  }

  if("EIC"%in%grab_what){
    checkProvidedMzPpm(mz, ppm)
    if(verbosity>1)last_time <- timeReport(last_time, text = "Extracting EIC...")
    if(!"MS1"%in%grab_what){
      init_dt <- grabMzxmlMS1(xml_data = xml_data, file_metadata = file_metadata,
                              rtrange = rtrange, prefilter = prefilter)
    } else {
      init_dt <- output_data$MS1
      if(!nrow(init_dt))stop("Something weird - can't find MS1 data to subset")
    }
    EIC_list <- lapply(unique(mz), function(mass){
      init_dt[mz%between%pmppm(mass = mass, ppm = ppm)]
    })
    output_data$EIC <- rbindlist(EIC_list)
  }

  if("EIC_MS2"%in%grab_what){
    checkProvidedMzPpm(mz, ppm)
    if(verbosity>1){
      last_time <- timeReport(last_time, text = "Extracting EIC MS2...")
      }
    if(!"MS2"%in%grab_what){
      init_dt <- grabMzxmlMS2(xml_data=xml_data, file_metadata = file_metadata,
                              rtrange=rtrange)
    } else {
      init_dt <- output_data$MS2
    }
    premz <- NULL #To prevent R CMD check "notes"
    EIC_MS2_list <- lapply(unique(mz), function(mass){
      init_dt[premz%between%pmppm(mass = mass, ppm = ppm)]
    })
    output_data$EIC_MS2 <- rbindlist(EIC_MS2_list)
  }

  if("metadata"%in%grab_what){
    if(verbosity>1){
      last_time <- timeReport(last_time, text = "Reading file metadata...")
    }
    output_data$metadata <- grabMzxmlMetadata(xml_data = xml_data)
  }

  if(verbosity>1){
    time_total <- round(difftime(Sys.time(), last_time), digits = 2)
    cat(time_total, units(time_total), "\n")
  }

  output_data
}



# Get mzXML specifics (functions of xml_data) ----

#' Helper function to extract mzXML file metadata
#'
#' @param xml_data mzXML data as parsed by xml2
#'
#' @return A list of values corresponding to various pieces of metadata
#' for each file
grabMzxmlMetadata <- function(xml_data){
  source_node <- xml2::xml_find_first(xml_data, xpath = "//d1:parentFile")
  if(length(source_node)>0){
    source_file <- basename(xml2::xml_attr(source_node, "fileName"))
  } else {
    source_file <- "None found"
  }

  inst_xpath <- "//d1:msInstrument/child::node()[starts-with(name(), 'ms')]"
  inst_nodes <- xml2::xml_find_all(xml_data, xpath = inst_xpath)
  if(length(inst_nodes)>0){
    inst_names <- xml2::xml_attr(inst_nodes, "category")
    inst_vals <- xml2::xml_attr(inst_nodes, "value")
    names(inst_vals) <- inst_names
  } else {
    inst_vals <- "None found"
    names(inst_vals) <- "Instrument data"
  }

  scan_nodes <- xml2::xml_find_all(xml_data, xpath = "//d1:scan")

  n_spectra <- length(scan_nodes)

  if(n_spectra>0){
    centroided <- as.integer(unique(xml2::xml_attr(scan_nodes, "centroided")))
    # 1 if centroided, 0 if not (i.e. profile mode)
    if (1 %in% centroided) {
      centroided <- TRUE
    } else {
      centroided <- FALSE
    }

    ms_levels <- paste0("MS", unique(xml2::xml_attr(scan_nodes, "msLevel")),
                        collapse = ", ")

    mz_lowest <- min(as.numeric(xml2::xml_attr(scan_nodes, "lowMz")))

    mz_highest <- max(as.numeric(xml2::xml_attr(scan_nodes, "highMz")))

    rt <- xml2::xml_attr(scan_nodes, "retentionTime")
    unit <- unique(gsub(".*[0-9]", "", rt))
    rt <- gsub("[^0-9.-]", "", rt)
    rt <- as.numeric(rt)
    if("S"%in%unit) rt <- rt/60

    rt_start <- min(rt)
    rt_end <- max(rt)

    polarities <- unique(xml2::xml_attr(scan_nodes, "polarity"))
    polarities[polarities %in% "+"] <- "positive"
    polarities[polarities %in% "-"] <- "negative"

  } else {
    centroided <- NA
    ms_levels <- NA_integer_
    mz_lowest <- NA_real_
    mz_highest <- NA_real_
    rt_start <- NA_real_
    rt_end <- NA_real_
    polarities <- NA_character_
  }

  metadata <- data.table(
    source_file=source_file,
    inst_data=list(inst_vals),
    n_spectra=n_spectra,
    ms_levels=ms_levels,
    mz_lowest=mz_lowest,
    mz_highest=mz_highest,
    rt_start=rt_start,
    rt_end=rt_end,
    centroided=centroided,
    polarity=polarities
  )
}


#' Helper function to extract mzXML file metadata
#'
#' @param xml_data mzXML data as parsed by xml2
#'
#' @return A list of values used by other parsing functions, currently
#' compression, precision, and endian encoding (endi_enc)
grabMzxmlEncodingData <- function(xml_data){
  peak_metadata <- xml2::xml_find_first(xml_data, '//d1:peaks')
  compr_type <- xml2::xml_attr(peak_metadata, "compressionType")
  compr <- switch(compr_type,
                  `zlib`="gzip",
                  `zlib compression`="gzip",
                  `no compression`="none",
                  `none`="none")
  if(is.null(compr))compr<-"none"

  enc_type <- xml2::xml_attr(peak_metadata, "precision")
  precision <- as.numeric(enc_type)/8

  byte_order <- xml2::xml_attr(peak_metadata, "byteOrder")
  endi_enc <- switch(byte_order, `network`="big")

  list(compression=compr, precision=precision, endi_enc=endi_enc)
}


#' Extract the MS1 data from an mzXML nodeset
#'
#' @param xml_data An `xml2` nodeset, usually created by applying `read_xml` to
#'   an mzXML file.
#' @param file_metadata Information about the file used to decode the binary
#'   arrays containing m/z and intensity information.
#' @param rtrange A vector of length 2 containing an upper and lower bound on
#'   retention times of interest. Providing a range here can speed up load times
#'   (although not enormously, as the entire file must still be read) and reduce
#'   the final object's size.
#' @param prefilter The lowest intensity value of interest, used to reduce file
#'   size (and especially useful for profile mode data with many 0 values)
#'
#' @return A `data.table` with columns for retention time (rt), m/z (mz),
#' and intensity (int).
grabMzxmlMS1 <- function(xml_data, file_metadata, rtrange, prefilter){
  ms1_xpath <- '//d1:scan[@msLevel="1" and @peaksCount>0]'
  ms1_nodes <- xml2::xml_find_all(xml_data, ms1_xpath)

  rt_vals <- grabMzxmlSpectraRt(ms1_nodes)
  if(!is.null(rtrange)){
    ms1_nodes <- ms1_nodes[rt_vals%between%rtrange]
    rt_vals <- rt_vals[rt_vals%between%rtrange]
  }

  mz_int_vals <- grabMzxmlSpectraMzInt(ms1_nodes, file_metadata)

  dt_data <- mapply(cbind, rt_vals, mz_int_vals, SIMPLIFY = FALSE)
  dt <- as.data.table(do.call(what=rbind, dt_data))
  names(dt) <- c("rt", "mz", "int")
  int <- NULL #To prevent R CMD check "notes" when using data.table syntax
  dt[int>prefilter]
}


#' Extract the MS2 data from an mzXML nodeset
#'
#' @param xml_data An `xml2` nodeset, usually created by applying `read_xml` to
#'   an mzXML file.
#' @param file_metadata Information about the file used to decode the binary
#'   arrays containing m/z and intensity information.
#' @param rtrange A vector of length 2 containing an upper and lower bound on
#'   retention times of interest. Providing a range here can speed up load times
#'   (although not enormously, as the entire file must still be read) and reduce
#'   the final object's size.
#'
#' @return A `data.table` with columns for retention time (rt),  precursor m/z (mz),
#' fragment m/z (fragmz), collision energy (voltage), and intensity (int).
grabMzxmlMS2 <- function(xml_data, file_metadata, rtrange){
  ms2_xpath <- '//d1:scan[@msLevel="2" and @peaksCount>0]'
  ms2_nodes <- xml2::xml_find_all(xml_data, ms2_xpath)
  if(!length(ms2_nodes)){
    return(data.table(rt=numeric(), premz=numeric(), fragmz=numeric(),
                      int=numeric(), voltage=integer()))
  }
  # Remove all nodes with zero peaks
  ms2_nodes <- ms2_nodes[as.numeric(xml2::xml_attr(ms2_nodes, "peaksCount"))>0]

  rt_vals <- grabMzxmlSpectraRt(ms2_nodes)
  if(!is.null(rtrange)){
    ms2_nodes <- ms2_nodes[rt_vals%between%rtrange]
    rt_vals <- rt_vals[rt_vals%between%rtrange]
  }

  premz_vals <- grabMzxmlSpectraPremz(ms2_nodes)
  voltage_vals <- grabMzxmlSpectraVoltage(ms2_nodes)
  mz_int_vals <- grabMzxmlSpectraMzInt(ms2_nodes, file_metadata)

  dt_data <- mapply(cbind, rt_vals, premz_vals, mz_int_vals,
                    voltage_vals, SIMPLIFY = FALSE)
  dt <- as.data.table(do.call(what=rbind, dt_data))
  names(dt) <- c("rt", "premz", "fragmz", "int", "voltage")
  dt$voltage <- as.integer(dt$voltage)
  dt
}


#' Grab the BPC or TIC from a file
#'
#' The base peak intensity and total ion current are actually written into the
#' mzXML files and aren't encoded, making retrieval of BPC and TIC information
#' blazingly fast if parsed correctly.
#'
#' @param xml_data An `xml2` nodeset, usually created by applying `read_xml` to
#'   an mzML file.
#' @param TIC Boolean. If TRUE, the TIC is extracted rather than the BPC.
#' @param rtrange A vector of length 2 containing an upper and lower bound on
#'   retention times of interest. Providing a range here can speed up load times
#'   (although not enormously, as the entire file must still be read) and reduce
#'   the final object's size.
#'
#' @return A `data.table` with columns for retention time (rt), and intensity (int).
grabMzxmlBPC <- function(xml_data, TIC=FALSE, rtrange){
  scan_nodes <- xml2::xml_find_all(
    xml_data, '//d1:scan[@msLevel="1"]'
  )
  rt_chrs <- xml2::xml_attr(scan_nodes, "retentionTime")
  rt_vals <- as.numeric(gsub(pattern = "PT|S", replacement = "", rt_chrs))
  if(any(rt_vals>150)){rt_vals <- rt_vals/60}

  int_attr <- ifelse(TIC, "totIonCurrent", "basePeakIntensity")
  int_vals <- as.numeric(xml2::xml_attr(scan_nodes, int_attr))
  if(!is.null(rtrange)){
    int_vals <- int_vals[rt_vals%between%rtrange]
    rt_vals <- rt_vals[rt_vals%between%rtrange]
  }


  return(data.table(rt=rt_vals, int=int_vals))
}


# Get spectrum things (functions of xml_nodes) ----

#' Extract the retention time from the spectra of an mzXML nodeset
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#' by the mass spectrometer, usually produced by applying `xml_find_all` to an
#' MS1 or MS2 nodeset.
#'
#' @return A numeric vector of retention times, one for each scan
#'

grabMzxmlSpectraRt <- function(xml_nodes){
  rt_attrs <- xml2::xml_attr(xml_nodes, "retentionTime")
  rt_unit <- unique(gsub(".*[0-9]", "", rt_attrs))
  rt_vals <- as.numeric(gsub("PT|S", "", rt_attrs))

  if ("S" %in% rt_unit) rt_vals <- rt_vals/60
  # if(any(rt_vals>150)){
  #   # Guess RT is in seconds if the run is more than 150 long
  #   # A 2.5 minute run is unheard of, and a 2.5 hour run is unheard of
  #   rt_vals <- rt_vals/60
  # }

  rt_vals
}


#' Extract the precursor mass from the spectra of an mzXML nodeset
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#' by the mass spectrometer, usually produced by applying `xml_find_all` to an
#' MS1 or MS2 nodeset.
#'
#' @return A numeric vector of precursor masses, one for each scan
#'

grabMzxmlSpectraPremz <- function(xml_nodes){
  premz_nodes <- xml2::xml_find_all(xml_nodes, xpath = "d1:precursorMz")
  as.numeric(xml2::xml_text(premz_nodes))
}


#' Extract the collison energies from the spectra of an mzXML nodeset
#'
#' Although the collision energy is typically fixed per file, it's equally
#' fast (afaik) to just grab them all individually here. Also, I'm worried about
#' these rumors of "ramped" collision energies
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#' by the mass spectrometer, usually produced by applying `xml_find_all` to an
#' MS1 or MS2 nodeset.
#'
#' @return A numeric vector of collision energies, one for each scan.
#'

grabMzxmlSpectraVoltage <- function(xml_nodes){
  as.numeric(xml2::xml_attr(xml_nodes, "collisionEnergy"))
}



#' Extract the mass-to-charge data from the spectra of an mzXML nodeset
#'
#' The mz and intensity information of mzXML files are encoded as a binary
#' array, sometimes compressed via gzip or zlib or numpress. This code finds all
#' the m/z-int binary arrays and converts them back to the original
#' measurements. See https://github.com/ProteoWizard/pwiz/issues/1301
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#'   by the mass spectrometer, usually produced by applying `xml_find_all` to an
#'   MS1 or MS2 nodeset.
#' @param file_metadata Information about the file used to decode the binary
#'   arrays containing m/z and intensity information. Here, the compression and
#'   mz precision information is relevant.
#'
#' @return A numeric vector of masses, many for each scan.
#'

grabMzxmlSpectraMzInt <- function(xml_nodes, file_metadata){
  all_peak_nodes <- xml2::xml_text(xml2::xml_find_all(xml_nodes, xpath = "d1:peaks"))
  vals <- lapply(all_peak_nodes, function(binary){
    if(!nchar(binary))return(matrix(ncol = 2, nrow = 0))
    decoded_binary <- base64enc::base64decode(binary)
    raw_binary <- as.raw(decoded_binary)
    decomp_binary <- memDecompress(raw_binary, type = file_metadata$compression)
    final_binary <- readBin(decomp_binary, what = "numeric",
                            n=length(decomp_binary)/file_metadata$precision,
                            size = file_metadata$precision,
                            endian = file_metadata$endi_enc)
    matrix(final_binary, ncol = 2, byrow = TRUE)
  })
}


# Other helper functions ----
shrinkRTrangemzXML <- function(xml_nodes, rtrange){
  rt_xpath <- 'd1:scanList/d1:scan/d1:cvParam[@name="scan start time"]'
  rt_nodes <- xml2::xml_find_all(xml_nodes, rt_xpath)
  rt_vals <- as.numeric(xml2::xml_attr(rt_nodes, "value"))
  if(any(rt_vals>150)){
    # Guess RT is in seconds if the run is more than 150 long
    # A 2.5 minute run is unheard of, and a 2.5 hour run is unheard of
    rt_vals <- rt_vals/60
  }
  xml_nodes[rt_vals%between%rtrange]
}
