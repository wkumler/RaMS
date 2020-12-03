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
#' @param verbose Boolean. If TRUE, R will print information about timing to the
#'   console as the file is read in.
#' @param mz A vector of the mass-to-charge ratio for compounds of interest.
#'   Only used when combined with `grab_what = "EIC"` (see above). Multiple
#'   masses can be provided.
#' @param ppm A single number corresponding to the mass accuracy (in parts per
#'   million) of the instrument on which the data was collected. Only used when
#'   combined with `grab_what = "EIC"` (see above).
#' @param rtrange Not supported for mzXML data. Only provided here so as to
#'   throw a friendly warning rather than an unexpected error.
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
#' sample_file <- system.file("extdata", "FK180310_Full1.mzXML.gz",
#'                            package = "RaMS")
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
#' sample_file <- system.file("extdata", "FK180310_DDApos100.mzML.gz",
#'                            package = "RaMS")
#' MS2_data <- grabMzmlData(sample_file, grab_what="MS2")
grabMzxmlData <- function(filename, grab_what, verbose=FALSE,
                          rtrange=NULL, mz=NULL, ppm=NULL){
  if(!is.null(rtrange)){
    message("\n`rtrange` argument not supported for mzXML files, ignoring")
  }
  if(verbose){
    cat(paste0("\nReading file ", basename(filename), "... "))
    last_time <- Sys.time()
  }
  xml_data <- xml2::read_xml(filename)

  checkFileType(xml_data, "mzXML")
  file_metadata <- grabMzxmlEncodingData(xml_data)

  output_data <- list()

  if("everything"%in%grab_what){
    if(length(setdiff(grab_what, "everything"))&&verbose){
      message(paste("Heads-up: grab_what = `everything` includes",
                    "MS1, MS2, BPC, and TIC data"))
      message("Ignoring additional grabs")
    }
    grab_what <- c("MS1", "MS2", "BPC", "TIC", "metadata")
  }

  if("MS1"%in%grab_what){
    if(verbose)last_time <- timeReport(last_time, text = "Reading MS1 data...")
    output_data$MS1 <- grabMzxmlMS1(xml_data = xml_data,
                                    file_metadata = file_metadata)
  }

  if("MS2"%in%grab_what){
    if(verbose)last_time <- timeReport(last_time, text = "Reading MS2 data...")
    output_data$MS2 <- grabMzxmlMS2(xml_data = xml_data,
                                    file_metadata = file_metadata)
  }

  if("BPC"%in%grab_what){
    if(verbose)last_time <- timeReport(last_time, text = "Reading BPC...")
    output_data$BPC <- grabMzxmlBPC(xml_data = xml_data)
  }

  if("TIC"%in%grab_what){
    if(verbose)last_time <- timeReport(last_time, text = "Reading TIC...")
    output_data$TIC <- grabMzxmlBPC(xml_data = xml_data, TIC = TRUE)
  }

  if("EIC"%in%grab_what){
    checkProvidedMzPpm(mz, ppm)
    if(verbose)last_time <- timeReport(last_time, text = "Extracting EIC...")
    if(!"MS1"%in%grab_what){
      init_dt <- grabMzxmlMS1(xml_data=xml_data, file_metadata = file_metadata)
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
    if(verbose){
      last_time <- timeReport(last_time, text = "Extracting EIC MS2...")
      }
    if(!"MS2"%in%grab_what){
      init_dt <- grabMzxmlMS2(xml_data=xml_data, file_metadata = file_metadata)
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
    if(verbose){
      last_time <- timeReport(last_time, text = "Reading file metadata...")
    }
    output_data$metadata <- grabMzxmlMetadata(xml_data = xml_data)
  }

  if(verbose){
    cat(Sys.time()-last_time, "s\n")
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

  mslevel_nodes <- xml2::xml_find_all(xml_data, xpath = "//d1:scan")
  if(length(mslevel_nodes)>0){
    mslevels <- paste0("MS", unique(xml_attr(mslevel_nodes, "msLevel")), collapse = ", ")
  } else {
    mslevels <- "None found"
  }


  metadata <- data.table(
    source_file=source_file,
    inst_data=list(inst_vals),
    mslevels=mslevels
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
#'
#' @return A `data.table` with columns for retention time (rt), m/z (mz),
#' and intensity (int).
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
grabMzxmlMS2 <- function(xml_data, file_metadata){
  ms2_xpath <- '//d1:scan[@msLevel="2"]'
  ms2_nodes <- xml2::xml_find_all(xml_data, ms2_xpath)
  if(!length(ms2_nodes)){
    return(data.table(rt=numeric(), premz=numeric(), fragmz=numeric(),
                      int=numeric(), voltage=integer()))
  }
  # Remove all nodes with zero peaks
  ms2_nodes <- ms2_nodes[as.numeric(xml2::xml_attr(ms2_nodes, "peaksCount"))>0]

  rt_vals <- grabMzxmlSpectraRt(ms2_nodes)
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
#'
#' @return A `data.table` with columns for retention time (rt), and intensity (int).
grabMzxmlBPC <- function(xml_data, TIC=FALSE){
  scan_nodes <- xml2::xml_find_all(xml_data, '//d1:scan[@msLevel="1"]')
  rt_chrs <- xml2::xml_attr(scan_nodes, "retentionTime")
  rt_vals <- as.numeric(gsub(pattern = "PT|S", replacement = "", rt_chrs))/60
  if(!any(rt_vals>1))stop("Are your mzXML files in minutes or seconds? (Should be sec)")

  int_attr <- ifelse(TIC, "totIonCurrent", "basePeakIntensity")
  int_vals <- as.numeric(xml2::xml_attr(scan_nodes, int_attr))

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
  rt_vals <- as.numeric(gsub("PT|S", "", rt_attrs))/60 # Convert to minutes
  if(!any(rt_vals>1))stop("Are your mzXML files in minutes or seconds? (Should be sec)")
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


