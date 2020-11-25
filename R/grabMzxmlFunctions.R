# Welcome to RaMS!

# grabMzxmlData ----

#' Get mass-spectrometry data from an mzXML file
#'
#' This function handles the mzXML side of things, reading in files that are
#' written in the mzXML format. Much of the code is similar to the mzXML format,
#' but the xpath handles are different and the mz/int array is encoded
#' simultaneously rather than as two separate entries.
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
#' @param verbose Boolean. If TRUE, R will print information about timing to the
#'   console as the file is read in.
#' @param mz A vector of the mass-to-charge ratio for compounds of interest.
#'   Only used when combined with `grab_what = "EIC"` (see above). Multiple
#'   masses can be provided.
#' @param ppm A single number corresponding to the mass accuracy (in parts per
#'   million) of the instrument on which the data was collected. Only used when
#'   combined with `grab_what = "EIC"` (see above).
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
grabMzxmlData <- function(filename, grab_what, verbose=FALSE,
                         mz=NULL, ppm=NULL){
  if(verbose){
    start_time <- Sys.time()
    last_time <- Sys.time()
    cat("\nReading file", basename(filename), "... ")
  }
  xml_data <- read_xml(filename)

  checkFileType(xml_data, "mzXML")
  file_metadata <- grabMzxmlEncodingData(xml_data)

  output_data <- list()

  if("everything"%in%grab_what){
    if(length(setdiff(grab_what, "everything"))&&verbose){
      message("Heads-up: grab_what = `everything` includes MS1, MS2, BPC, and TIC data")
      message("Ignoring additional grabs")
    }
    grab_what <- c("MS1", "MS2", "BPC", "TIC")
  }

  if("MS1"%in%grab_what){
    if(verbose){
      cat(Sys.time()-last_time, "s\n")
      last_time <- Sys.time()
      cat("Reading MS1 data... ")
    }
    output_data$MS1 <- grabMzxmlMS1(xml_data, file_metadata)
  }

  if("MS2"%in%grab_what){
    if(verbose){
      cat(Sys.time()-last_time, "s\n")
      last_time <- Sys.time()
      cat("Reading MS2 data... ")
    }
    output_data$MS2 <- grabMzxmlMS2(xml_data, file_metadata)
  }

  if("BPC"%in%grab_what){
    if(verbose){
      cat(Sys.time()-last_time, "s\n")
      last_time <- Sys.time()
      cat("Reading BPC... ")
    }
    output_data$BPC <- grabMzxmlBPC(xml_data)
  }

  if("TIC"%in%grab_what){
    if(verbose){
      cat(Sys.time()-last_time, "s\n")
      last_time <- Sys.time()
      cat("Reading TIC... ")
    }
    output_data$TIC <- grabMzxmlBPC(xml_data, TIC = TRUE)
  }

  # if("EIC"%in%grab_what){
  #   checkProvidedMzPpm(mz, ppm)
  #
  # }
  if(verbose){
    cat(Sys.time()-last_time, "s\n")
    cat("Total time:", Sys.time()-start_time, "\n")
  }

  output_data
}

# Get mzXML specifics (functions of xml_data) ----

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
#'
#' @return A `data.table` with columns for retention time (rt), m/z (mz), and intensity (int).
#'
#' @examples
grabMzxmlMS1 <- function(xml_data, file_metadata){
  ms1_xpath <- '//d1:scan[@msLevel="1"]'
  ms1_nodes <- xml2::xml_find_all(xml_data, ms1_xpath)

  rt_vals <- grabMzxmlSpectraRt(ms1_nodes)
  mz_int_vals <- grabMzxmlSpectraMzInt(ms1_nodes, file_metadata)

  dt_data <- mapply(cbind, rt_vals, mz_int_vals, SIMPLIFY = FALSE)
  dt <- as.data.table(do.call(what=rbind, dt_data))
  names(dt) <- c("rt", "mz", "int")
  dt
}


#' Extract the MS2 data from an mzXML nodeset
#'
#' @param xml_data An `xml2` nodeset, usually created by applying `read_xml` to
#'   an mzXML file.
#' @param file_metadata Information about the file used to decode the binary
#'   arrays containing m/z and intensity information.
#'
#' @return A `data.table` with columns for retention time (rt),  precursor m/z (mz),
#' fragment m/z (fragmz), collision energy (voltage), and intensity (int).
#'
#' @examples
grabMzxmlMS2 <- function(xml_data, file_metadata){
  ms2_xpath <- '//d1:scan[@msLevel="2"]'
  ms2_nodes <- xml2::xml_find_all(xml_data, ms2_xpath)
  if(!length(ms2_nodes)){
    return(data.table(rt=numeric(), premz=numeric(), fragmz=numeric(),
                      int=numeric(), voltages=numeric()))
  }

  rt_vals <- grabMzxmlSpectraRt(ms2_nodes)
  premz_vals <- grabMzxmlSpectraPremz(ms2_nodes)
  voltage_vals <- grabMzxmlSpectraVoltage(ms2_nodes)
  mz_int_vals <- grabMzxmlSpectraMzInt(ms2_nodes, file_metadata)

  dt_data <- mapply(cbind, rt_vals, premz_vals, mz_int_vals,
                    voltage_vals, SIMPLIFY = FALSE)
  dt <- as.data.table(do.call(what=rbind, dt_data))
  names(dt) <- c("rt", "premz", "fragmz", "int", "voltages")
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
#'
#' @return A `data.table` with columns for retention time (rt), and intensity (int).
#'
#' @examples
grabMzxmlBPC <- function(xml_data, TIC=FALSE){
  scan_nodes <- xml2::xml_find_all(xml_data, '//d1:scan[@msLevel="1"]')
  rt_chrs <- xml2::xml_attr(scan_nodes, "retentionTime")
  rt_vals <- as.numeric(gsub(pattern = "PT|S", replacement = "", rt_chrs))/60

  int_attr <- ifelse(TIC, "totIonCurrent", "basePeakIntensity")
  int_vals <- as.numeric(xml2::xml_attr(scan_nodes, int_attr))

  return(data.frame(rt=rt_vals, int=int_vals))
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
#' @examples
grabMzxmlSpectraRt <- function(xml_nodes){
  rt_attrs <- xml2::xml_attr(xml_nodes, "retentionTime")
  as.numeric(gsub("PT|S", "", rt_attrs))
}


#' Extract the precursor mass from the spectra of an mzXML nodeset
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#' by the mass spectrometer, usually produced by applying `xml_find_all` to an
#' MS1 or MS2 nodeset.
#'
#' @return A numeric vector of precursor masses, one for each scan
#'
#' @examples
grabMzxmlSpectraPremz <- function(xml_nodes){
  premz_nodes <- xml_find_all(xml_nodes, xpath = "d1:precursorMz")
  as.numeric(xml_text(premz_nodes))
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
#' @examples
grabMzxmlSpectraVoltage <- function(xml_nodes){
  filterline_data <- xml_attr(xml_nodes, "filterLine")
  as.numeric(gsub(".*cid| .*", "", filterline_data))
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
#' @examples
grabMzxmlSpectraMzInt <- function(xml_nodes, file_metadata){
  all_peak_nodes <- xml_text(xml_find_all(xml_nodes, xpath = "d1:peaks"))
  vals <- lapply(all_peak_nodes, function(binary){
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


