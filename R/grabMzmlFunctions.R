# Welcome to RaMS!

# grabMzmlData ----

#' Get mass-spectrometry data from an mzML file
#'
#' This function handles the mzML side of things, reading in files that are
#' written in the mzML format. Much of the code is similar to the mzXML format,
#' but the xpath handles are different and the mz/int array is encoded as two
#' separate entries rather than simultaneously. This function has been exposed
#' to the user in case per-file optimization (such as peakpicking or additional
#' filtering) is desired before the full data object is returned.
#'
#' @param filename A single filename to read into R's memory. Both absolute and
#'   relative paths are acceptable.
#' @param grab_what What data should be read from the file? Options include
#'   "MS1" for data only from the first spectrometer, "MS2" and "MS3" for
#'   fragmentation data, "BPC" for rapid access to the base peak chromatogram,
#'   "TIC" for rapid access to the total ion chromatogram, "DAD" for DAD (UV)
#'   data, and "chroms" for precompiled chromatogram data (especially useful for
#'   MRM but often contains BPC/TIC in other files). Metadata can be accessed
#'   with "metadata", which provides information about the instrument and time
#'   the file was run. These options can be combined (i.e. `grab_data=c("MS1",
#'   "MS2", "BPC")`) or this argument can be set to "everything" to extract all
#'   of the above. Options "EIC", "EIC_MS2", and "EIC_MS3" are useful when
#'   working with files whose total size exceeds working memory - it first
#'   extracts all relevant MS1/2/3 data, respectively, then discards data
#'   outside of the mass range(s) calculated from the provided mz and ppm. The
#'   default, "everything", includes all MS1, MS2, BPC, TIC, and metadata.
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
#' @param rtrange A vector of length 2 containing an upper and lower bound on
#'   retention times of interest. Providing a range here can speed up load times
#'   (although not enormously, as the entire file must still be read) and reduce
#'   the final object's size.
#' @param prefilter A single number corresponding to the minimum intensity of
#'   interest in the MS1 data. Data points with intensities below this threshold
#'   will be silently dropped, which can dramatically reduce the size of the
#'   final object. Currently only works with MS1 data, but could be expanded
#'   easily to handle more.
#' @param incl_polarity Toggle this option to TRUE for mixed-polarity files. An
#'   additional column will be added corresponding to the polarity of the scan,
#'   with either a 1 or a -1 corresponding to positive and negative mode,
#'   respectively.
#'
#' @return A list of `data.table`s, each named after the arguments requested in
#'   grab_what. E.g. $MS1 contains MS1 information, $MS2 contains fragmentation
#'   info, etc. MS1 data has four columns: retention time (rt), mass-to-charge
#'   (mz), intensity (int), and filename. MS2 data has six: retention time (rt),
#'   precursor m/z (premz), fragment m/z (fragmz), fragment intensity (int),
#'   collision energy (voltage), and filename. MS3 has an additional column to
#'   MS2 (prepremz) which has the original MS1 scan's m/z ratio. Data requested
#'   that does not exist in the provided files (such as MS2 data requested from
#'   MS1-only files) will return an empty (length zero) data.table. The
#'   data.tables extracted from each of the individual files are collected into
#'   one large table using data.table's `rbindlist`. $metadata is a little
#'   weirder because the metadata doesn't fit neatly into a tidy format but
#'   things are hopefully named helpfully. $chroms was added in v1.3 and
#'   contains 7 columns: chromatogram type (usually TIC, BPC or SRM info),
#'   chromatogram index, target mz, product mz, retention time (rt), and
#'   intensity (int). $DAD was also added in v1.3 and contains has three
#'   columns: retention time (rt), wavelength (lambda),and intensity (int).
#'
#' @export
#'
#' @examples
#' \dontshow{data.table::setDTthreads(2)}
#' sample_file <- system.file("extdata", "LB12HL_AB.mzML.gz", package = "RaMS")
#' file_data <- grabMzmlData(sample_file, grab_what="MS1")
#' \dontrun{
#' # Extract MS1 data and a base peak chromatogram
#' file_data <- grabMzmlData(sample_file, grab_what=c("MS1", "BPC"))
#' # Extract data from a retention time subset
#' file_data <- grabMzmlData(sample_file, grab_what=c("MS1", "BPC"),
#'                           rtrange=c(5, 7))
#' # Extract EIC for a specific mass
#' file_data <- grabMzmlData(sample_file, grab_what="EIC", mz=118.0865, ppm=5)
#' # Extract EIC for several masses simultaneously
#' file_data <- grabMzmlData(sample_file, grab_what="EIC", ppm=5,
#'                           mz=c(118.0865, 146.118104, 189.123918))
#'
#' # Extract MS2 data
#' sample_file <- system.file("extdata", "S30657.mzML.gz", package = "RaMS")
#' MS2_data <- grabMzmlData(sample_file, grab_what="MS2")
#' }
grabMzmlData <- function(filename, grab_what, verbosity=0, incl_polarity=FALSE,
                         mz=NULL, ppm=NULL, rtrange=NULL, prefilter=-1){
  if(verbosity>1){
    cat(paste0("\nReading file ", basename(filename), "... "))
    last_time <- Sys.time()
  }
  xml_data <- xml2::read_xml(filename)

  checkNamespace(xml_data)
  checkFileType(xml_data, "mzML")
  rtrange <- checkRTrange(rtrange)
  prefilter <- checkProvidedPrefilter(prefilter)

  output_data <- list()

  if("everything"%in%grab_what){
    extra_grabs <- setdiff(grab_what, "everything")
    if(any(c("MS1", "MS2", "BPC", "TIC", "metadata")%in%extra_grabs)&&verbosity>0){
      message(paste("Heads-up: grab_what = `everything` includes",
                    "MS1, MS2, BPC, TIC, and meta data"))
      message("Ignoring duplicate specification")
    }
    grab_what <- unique(c("MS1", "MS2", "BPC", "TIC", "metadata", extra_grabs))
  }

  if(any(c("MS1", "MS2", "MS3", "DAD", "EIC", "EIC_MS2", "EIC_MS3", "chroms")%in%grab_what)){
    file_metadata <- grabMzmlEncodingData(xml_data)
  }

  if("MS1"%in%grab_what){
    if(verbosity>1)last_time <- timeReport(last_time, text = "Reading MS1 data...")
    output_data$MS1 <- grabMzmlMS1(xml_data = xml_data, rtrange = rtrange,
                                   file_metadata = file_metadata,
                                   incl_polarity = incl_polarity,
                                   prefilter = prefilter)
  }

  if("MS2"%in%grab_what){
    if(verbosity>1)last_time <- timeReport(last_time, text = "Reading MS2 data...")
    output_data$MS2 <- grabMzmlMS2(xml_data = xml_data, rtrange = rtrange,
                                   incl_polarity = incl_polarity,
                                   file_metadata = file_metadata)
  }

  if("MS3"%in%grab_what){
    if(verbosity>1)last_time <- timeReport(last_time, text = "Reading MS3 data...")
    output_data$MS3 <- grabMzmlMS3(xml_data = xml_data, rtrange = rtrange,
                                   incl_polarity = incl_polarity,
                                   file_metadata = file_metadata)
  }

  if("DAD"%in%grab_what){
    if(verbosity>1)last_time <- timeReport(last_time, text = "Reading DAD data...")
    output_data$DAD <- grabMzmlDAD(xml_data = xml_data, rtrange = rtrange,
                                   file_metadata = file_metadata)
  }

  if("BPC"%in%grab_what){
    if(verbosity>1)last_time <- timeReport(last_time, text = "Reading BPC...")
    output_data$BPC <- grabMzmlBPC(xml_data = xml_data, rtrange = rtrange,
                                   incl_polarity = incl_polarity)
  }

  if("TIC"%in%grab_what){
    if(verbosity>1)last_time <- timeReport(last_time, text = "Reading TIC...")
    output_data$TIC <- grabMzmlBPC(xml_data = xml_data, rtrange = rtrange,
                                   incl_polarity = incl_polarity, TIC = TRUE)
  }

  if("EIC"%in%grab_what){
    checkProvidedMzPpm(mz, ppm)
    if(verbosity>1){
      last_time <- timeReport(last_time, text = "Extracting EIC...")
    }
    if(!"MS1"%in%grab_what){
      init_dt <- grabMzmlMS1(xml_data = xml_data, rtrange = rtrange,
                             file_metadata = file_metadata,
                             incl_polarity = incl_polarity,
                             prefilter = prefilter)
    } else {
      init_dt <- output_data$MS1
      if(!nrow(init_dt))stop("Something weird - can't find MS1 data to subset")
    }
    EIC_list <- lapply(unique(mz), function(mass){
      init_dt[mz%between%pmppm(mass = mass, ppm = ppm)]
    })
    output_data$EIC <- unique(rbindlist(EIC_list))
  }

  if("EIC_MS2"%in%grab_what){
    checkProvidedMzPpm(mz, ppm)
    if(verbosity>1){
      last_time <- timeReport(last_time, text = "Extracting EIC MS2...")
    }
    if(!"MS2"%in%grab_what){
      init_dt <- grabMzmlMS2(xml_data = xml_data, rtrange = rtrange,
                             file_metadata = file_metadata,
                             incl_polarity = incl_polarity)
    } else {
      init_dt <- output_data$MS2
    }
    premz <- NULL #To prevent R CMD check "notes"
    EIC_MS2_list <- lapply(unique(mz), function(mass){
      init_dt[premz%between%pmppm(mass = mass, ppm = ppm)]
    })
    output_data$EIC_MS2 <- unique(rbindlist(EIC_MS2_list))
  }

  if("EIC_MS3"%in%grab_what){
    checkProvidedMzPpm(mz, ppm)
    if(verbosity>1){
      last_time <- timeReport(last_time, text = "Extracting EIC MS3...")
    }
    if(!"MS3"%in%grab_what){
      init_dt <- grabMzmlMS3(xml_data = xml_data, rtrange = rtrange,
                             file_metadata = file_metadata,
                             incl_polarity = incl_polarity)
    } else {
      init_dt <- output_data$MS3
    }
    prepremz <- NULL #To prevent R CMD check "notes"
    EIC_MS3_list <- lapply(unique(mz), function(mass){
      init_dt[prepremz%between%pmppm(mass = mass, ppm = ppm)]
    })
    output_data$EIC_MS3 <- unique(rbindlist(EIC_MS3_list))
  }

  if("chroms"%in%grab_what){
    if(verbosity>1){
      last_time <- timeReport(last_time, text = "Reading chromatograms...")
    }
    output_data$chroms <- grabMzmlChroms(xml_data = xml_data, file_metadata = file_metadata)
  }

  if("metadata"%in%grab_what){
    if(verbosity>1){
      last_time <- timeReport(last_time, text = "Reading file metadata...")
    }
    output_data$metadata <- grabMzmlMetadata(xml_data = xml_data)
  }

  if(verbosity>1){
    time_total <- round(difftime(Sys.time(), last_time), digits = 2)
    cat(time_total, units(time_total), "\n")
  }

  output_data
}

# Get mzML specifics (functions of xml_data) ----

#' Helper function to extract mzML file metadata
#'
#' @param xml_data mzML data as parsed by xml2
#'
#' @return A list of values corresponding to various pieces of metadata
#' for each file
grabMzmlMetadata <- function(xml_data){
  source_node <- xml2::xml_find_first(xml_data, xpath = "//d1:sourceFile")
  if(length(source_node)>0){
    source_file <- xml2::xml_attr(source_node, "name")
  } else {
    source_file <- "None found"
  }

  inst_xpath <- "//d1:referenceableParamGroup/d1:cvParam"
  inst_nodes <- xml2::xml_find_first(xml_data, xpath = inst_xpath)
  if(length(inst_nodes)>0){
    inst_val <- xml2::xml_attr(inst_nodes, "name")
  } else {
    inst_val <- "None found"
  }

  config_xpath <- "//d1:componentList/child::node()"
  config_nodes <- xml2::xml_find_all(xml_data, xpath = config_xpath)
  if(length(inst_nodes)>0){
    config_types <- xml2::xml_name(config_nodes)
    config_order <- xml2::xml_attr(config_nodes, "order")
    config_name_nodes <- xml2::xml_find_first(config_nodes, "d1:cvParam")
    config_names <- xml2::xml_attr(config_name_nodes, "name")
  } else {
    config_types <- "None found"
    config_order <- "None found"
    config_names <- "None found"
  }


  time_node <- xml2::xml_find_first(xml_data, xpath = "//d1:run")
  time_val <- xml2::xml_attr(time_node, "startTimeStamp")
  if(!is.na(time_val)){
    time_stamp <- as.POSIXct(strptime(time_val, "%Y-%m-%dT%H:%M:%SZ"))
  } else {
    time_stamp <- as.POSIXct(NA)
  }

  mslevel_xpath <- '//d1:spectrum/d1:cvParam[@name="ms level"]'
  mslevel_nodes <- xml2::xml_find_all(xml_data, xpath = mslevel_xpath)
  if(length(mslevel_nodes)>0){
    ms_levels <- paste0("MS", unique(xml2::xml_attr(mslevel_nodes, "value")),
                      collapse = ", ")
  } else {
    ms_levels <- "None found"
  }

  mzlow_xpath <- '//d1:spectrum/d1:cvParam[@name="lowest observed m/z"]'
  mzlow_nodes <- xml2::xml_find_all(xml_data, xpath = mzlow_xpath)
  if(length(mzlow_nodes)>0){
    mz_lowest <- min(as.numeric(xml2::xml_attr(mzlow_nodes, "value")))
  } else {
    mz_lowest <- NA_real_
  }

  mzhigh_xpath <- '//d1:spectrum/d1:cvParam[@name="highest observed m/z"]'
  mzhigh_nodes <- xml2::xml_find_all(xml_data, xpath = mzhigh_xpath)
  if(length(mzhigh_nodes)>0){
    mz_highest <- max(as.numeric(xml2::xml_attr(mzhigh_nodes, "value")))
  } else {
    mz_highest <- NA_real_
  }

  rt_xpath <- '//d1:spectrum/d1:scanList/d1:scan/d1:cvParam[@name="scan start time"]'
  rt_nodes <- xml2::xml_find_all(xml_data, xpath = rt_xpath)
  rt_unit <- unique(xml_attr(rt_nodes, "unitName"))
  rt <- as.numeric(xml2::xml_attr(rt_nodes, "value"))

  if (!"minute" %in% rt_unit) rt=rt/60

  if(length(rt) > 0){
    rt_start <- min(rt)
    rt_end <- max(rt)
  } else {
    rt_start <- NA_real_
    rt_end <- NA_real_
  }

  centroided_xpath <- '//d1:spectrum/d1:cvParam[@accession="MS:1000127"]'
  centroided_nodes <- xml2::xml_find_all(xml_data, xpath = centroided_xpath)
  if (length(centroided_nodes) > 0) {
    centroided <- TRUE
  } else {
    profile_xpath <- '//d1:spectrum/d1:cvParam[@accession="MS:1000128"]'
    profile_nodes <- xml2::xml_find_all(xml_data, xpath = profile_xpath)
    if (length(profile_nodes) > 0) {
      centroided <- FALSE
    } else {
      centroided <- NA
    }
  }

  polarity_pos <- '//d1:spectrum/d1:cvParam[@accession="MS:1000130"]'
  polarity_pos <- xml_find_all(xml_data, polarity_pos)

  polarity_neg <- '//d1:spectrum/d1:cvParam[@accession="MS:1000129"]'
  polarity_neg <- xml_find_all(xml_data, polarity_neg)


  if(length(polarity_pos)>0|length(polarity_neg)>0) {
    polarities <- c(
      unique(gsub(" scan", "", xml_attr(polarity_pos, "name"))),
      unique(gsub(" scan", "", xml_attr(polarity_neg, "name")))
    )
  } else {
    polarities <- NA_character_
  }

  lambda_high_xpath <- '//d1:spectrum/d1:cvParam[@name="highest observed wavelength"]'
  lambda_high_nodes <- xml2::xml_find_all(xml_data, xpath = lambda_high_xpath)
  if(length(lambda_high_nodes)>0){
    lambda_highest <- max(as.numeric(xml2::xml_attr(lambda_high_nodes, "value")))
  } else {
    lambda_highest <- NA_real_
  }

  lambda_low_xpath <- '//d1:spectrum/d1:cvParam[@name="lowest observed wavelength"]'
  lambda_low_nodes <- xml2::xml_find_all(xml_data, xpath = lambda_low_xpath)
  if(length(lambda_low_nodes)>0){
    lambda_lowest <- max(as.numeric(xml2::xml_attr(lambda_low_nodes, "value")))
  } else {
    lambda_lowest <- NA_real_
  }

  n_spectra <- length(rt_nodes)

  chrom_xpath <- '//d1:chromatogram'
  chrom_nodes <- xml2::xml_find_all(xml_data, chrom_xpath)
  n_chromatograms <- length(chrom_nodes)

  metadata <- data.table(
    source_file=source_file,
    inst_data=inst_val,
    config_data=list(data.frame(
      order=config_order,
      type=config_types,
      name=config_names
    )),
    timestamp = time_stamp,
    n_spectra=n_spectra,
    n_chromatograms=n_chromatograms,
    ms_levels=ms_levels,
    mz_lowest=mz_lowest,
    mz_highest=mz_highest,
    lambda_lowest=lambda_lowest,
    lambda_highest=lambda_highest,
    rt_start=rt_start,
    rt_end=rt_end,
    centroided=centroided,
    polarity=polarities
  )
}

#' Helper function to extract mzML file encoding data
#'
#' @param xml_data mzML data as parsed by xml2
#'
#' @return A list of values used by other parsing functions, currently
#' compression, mz_precision, int_precision
grabMzmlEncodingData <- function(xml_data){
  init_xpath <- "//*[self::d1:spectrum or self::d1:chromatogram]"
  init_node <- xml2::xml_find_first(xml_data, xpath = init_xpath)
  if(length(init_node)==0){
    stop(paste("Unable to find a spectrum or chromatogram node from",
               "which to extract metadata"))
  }
  compr_xpath <- paste0('//d1:cvParam[@accession="MS:1000574"]|',
                        '//d1:cvParam[@accession="MS:1000576"]')
  compr_node <- xml2::xml_find_first(init_node, compr_xpath)
  compr_type <- xml2::xml_attr(compr_node, "name")
  compr <- switch(compr_type,
                  `zlib`="gzip",
                  `zlib compression`="gzip",
                  `no compression`="none",
                  `none`="none")

  mz_precision_xpath <- '//d1:cvParam[@accession="MS:1000523"]'
  mz_bit_node <- xml2::xml_find_first(init_node, mz_precision_xpath)
  mz_bit_type <- xml2::xml_attr(mz_bit_node, "name")
  mz_precision <- sub(mz_bit_type, pattern = "-bit float", replacement = "")
  mz_precision <- as.numeric(mz_precision)/8

  int_bit_xpath <- '//d1:cvParam[@accession="MS:1000521"]'
  int_bit_node <- xml2::xml_find_first(init_node, int_bit_xpath)
  int_bit_type <- xml2::xml_attr(int_bit_node, "name")
  int_precision <- sub(int_bit_type, pattern = "-bit float", replacement = "")
  int_precision <- as.numeric(int_precision)/8

  if(is.na(int_precision))int_precision <- mz_precision
  if(is.na(mz_precision))mz_precision <- int_precision

  list(compression=compr, mz_precision=mz_precision,
       int_precision=int_precision, endi_enc="little")
}


#' Extract the MS1 data from an mzML nodeset
#'
#' @param xml_data An `xml2` nodeset, usually created by applying `read_xml` to
#'   an mzML file.
#' @param rtrange A vector of length 2 containing an upper and lower bound on
#'   retention times of interest. Providing a range here can speed up load times
#'   (although not enormously, as the entire file must still be read) and reduce
#'   the final object's size.
#' @param file_metadata Information about the file used to decode the binary
#'   arrays containing m/z and intensity information.
#' @param prefilter The lowest intensity value of interest, used to reduce file
#'   size (and especially useful for profile mode data with many 0 values)
#'
#' @return A `data.table` with columns for retention time (rt), m/z (mz), and
#'   intensity (int).
grabMzmlMS1 <- function(xml_data, rtrange, file_metadata, prefilter, incl_polarity){
  ms1_xpath <- '//d1:spectrum[d1:cvParam[@name="ms level" and @value="1"]]'
  ms1_nodes <- xml2::xml_find_all(xml_data, ms1_xpath)
  if(!length(ms1_nodes)){
    if(incl_polarity){
      return(data.table(rt=numeric(), mz=numeric(), int=numeric(), polarity=numeric()))
    }
    return(data.table(rt=numeric(), mz=numeric(), int=numeric()))
  }

  rt_vals <- grabSpectraRt(ms1_nodes)
  if(!is.null(rtrange)){
    ms1_nodes <- ms1_nodes[rt_vals%between%rtrange]
    rt_vals <- rt_vals[rt_vals%between%rtrange]
  }

  mz_vals <- grabSpectraMz(ms1_nodes, file_metadata)
  int_vals <- grabSpectraInt(ms1_nodes, file_metadata)
  if(incl_polarity){
    pol_vals <- grabSpectraPolarity(ms1_nodes)
  }

  int <- NULL #To prevent R CMD check "notes"  when using data.table syntax
  if(incl_polarity){
    all_data <- data.table(rt=rep(rt_vals, sapply(mz_vals, length)),
                           mz=unlist(mz_vals), int=as.numeric(unlist(int_vals)),
                           polarity=rep(pol_vals, sapply(mz_vals, length)))
  } else {
    all_data <- data.table(rt=rep(rt_vals, sapply(mz_vals, length)),
                           mz=unlist(mz_vals), int=as.numeric(unlist(int_vals)))
  }
  all_data[int>prefilter]
}


#' Extract the MS2 data from an mzML nodeset
#'
#' @param xml_data An `xml2` nodeset, usually created by applying `read_xml` to
#'   an mzML file.
#' @param rtrange A vector of length 2 containing an upper and lower bound on
#'   retention times of interest. Providing a range here can speed up load times
#'   (although not enormously, as the entire file must still be read) and reduce
#'   the final object's size.
#' @param file_metadata Information about the file used to decode the binary
#'   arrays containing m/z and intensity information.
#'
#' @return A `data.table` with columns for retention time (rt),  precursor m/z
#'   (premz), fragment m/z (fragmz), collision energy (voltage), and intensity
#'   (int).
grabMzmlMS2 <- function(xml_data, rtrange, file_metadata, incl_polarity){
  ms2_xpath <- '//d1:spectrum[d1:cvParam[@name="ms level" and @value="2"]]'

  ms2_nodes <- xml2::xml_find_all(xml_data, ms2_xpath)
  if(!length(ms2_nodes)){
    return(data.table(rt=numeric(), premz=numeric(), fragmz=numeric(),
                      int=numeric(), voltage=integer()))
  }

  rt_vals <- grabSpectraRt(ms2_nodes)
  if(!is.null(rtrange)){
    ms2_nodes <- ms2_nodes[rt_vals%between%rtrange]
    rt_vals <- rt_vals[rt_vals%between%rtrange]
  }

  premz_vals <- grabSpectraPremz(ms2_nodes)
  voltage <- grabSpectraVoltage(ms2_nodes)
  if(all(is.na(voltage))){voltage <- rep(NA_integer_, length(premz_vals))}
  if(incl_polarity){pol_vals <- grabSpectraPolarity(ms2_nodes)}
  mz_vals <- grabSpectraMz(ms2_nodes, file_metadata)
  int_vals <- grabSpectraInt(ms2_nodes, file_metadata)

  if(length(premz_vals)==0 & length(voltage)==0 &
     length(mz_vals)==0 & length(int_vals)==0){
    if(incl_polarity){
      return(data.table(rt=numeric(), premz=numeric(), fragmz=numeric(),
                        int=numeric(), voltage=integer(), polarity=numeric()))
    } else {
      return(data.table(rt=numeric(), premz=numeric(), fragmz=numeric(),
                        int=numeric(), voltage=integer()))
    }
  }

  if(incl_polarity){
    return(
      data.table(rt=rep(rt_vals, sapply(mz_vals, length)),
                 premz=rep(premz_vals, sapply(mz_vals, length)),
                 fragmz=unlist(mz_vals), int=as.numeric(unlist(int_vals)),
                 voltage=rep(voltage, sapply(mz_vals, length)),
                 polarity=rep(pol_vals, sapply(mz_vals, length)))
    )
  }
  data.table(rt=rep(rt_vals, sapply(mz_vals, length)),
             premz=rep(premz_vals, sapply(mz_vals, length)),
             fragmz=unlist(mz_vals), int=as.numeric(unlist(int_vals)),
             voltage=rep(voltage, sapply(mz_vals, length)))
}


#' Extract the MS3 data from an mzML nodeset
#'
#' @param xml_data An `xml2` nodeset, usually created by applying `read_xml` to
#'   an mzML file.
#' @param rtrange A vector of length 2 containing an upper and lower bound on
#'   retention times of interest. Providing a range here can speed up load times
#'   (although not enormously, as the entire file must still be read) and reduce
#'   the final object's size.
#' @param file_metadata Information about the file used to decode the binary
#'   arrays containing m/z and intensity information.
#'
#' @return A `data.table` with columns for retention time (rt),
#' MS1 precursor m/z (prepremz), MS2 precursor m/z (premz),
#' fragment m/z (fragmz), collision energy (voltage), and intensity (int).
grabMzmlMS3 <- function(xml_data, rtrange, file_metadata, incl_polarity){
  ms3_xpath <- '//d1:spectrum[d1:cvParam[@name="ms level" and @value="3"]]'

  ms3_nodes <- xml2::xml_find_all(xml_data, ms3_xpath)
  if(!length(ms3_nodes)){
    return(data.table(rt=numeric(), prepremz=numeric(), premz=numeric(),
                      fragmz=numeric(), int=numeric(), voltage=integer()))
  }

  rt_vals <- grabSpectraRt(ms3_nodes)
  if(!is.null(rtrange)){
    ms3_nodes <- ms3_nodes[rt_vals%between%rtrange]
    rt_vals <- rt_vals[rt_vals%between%rtrange]
  }

  premz_vals <- matrix(grabSpectraPremz(ms3_nodes), ncol=2, byrow=TRUE)
  voltage <- matrix(grabSpectraVoltage(ms3_nodes), ncol=2, byrow = TRUE)
  if(all(is.na(voltage))){
    voltage <- matrix(rep(NA_integer_, nrow(premz_vals)*2), ncol=2, byrow=TRUE)
  }
  if(incl_polarity){pol_vals <- grabSpectraPolarity(ms3_nodes)}
  mz_vals <- grabSpectraMz(ms3_nodes, file_metadata)
  int_vals <- grabSpectraInt(ms3_nodes, file_metadata)

  if(length(premz_vals)==0 & length(voltage)==0 &
     length(mz_vals)==0 & length(int_vals)==0){
    if(incl_polarity){
      return(data.table(rt=numeric(), prepremz=numeric(), premz=numeric(),
                        fragmz=numeric(), int=numeric(), voltage=integer(),
                        polarity=numeric()))
    } else {
      return(data.table(rt=numeric(), prepremz=numeric(), premz=numeric(),
                        fragmz=numeric(), int=numeric(), voltage=integer()))
    }
  }

  if(incl_polarity){
    return(
      data.table(rt=rep(rt_vals, lengths(mz_vals)),
                 prepremz=rep(premz_vals[,2], lengths(mz_vals)),
                 premz=rep(premz_vals[,1], lengths(mz_vals)),
                 fragmz=unlist(mz_vals), int=as.numeric(unlist(int_vals)),
                 voltage=rep(voltage[,2], lengths(mz_vals)),
                 polarity=rep(pol_vals, lengths(mz_vals)))
    )
  }
  data.table(rt=rep(rt_vals, lengths(mz_vals)),
             prepremz=rep(premz_vals[,2], lengths(mz_vals)),
             premz=rep(premz_vals[,1], lengths(mz_vals)),
             fragmz=unlist(mz_vals), int=as.numeric(unlist(int_vals)),
             voltage=rep(voltage[,2], lengths(mz_vals)))
}


#' Grab the BPC or TIC from a file
#'
#' The base peak intensity and total ion current are actually written into the
#' mzML files and aren't encoded, making retrieval of BPC and TIC information
#' blazingly fast if parsed correctly.
#'
#' @param xml_data An `xml2` nodeset, usually created by applying `read_xml` to
#'   an mzML file.
#' @param rtrange A vector of length 2 containing an upper and lower bound on
#'   retention times of interest. Providing a range here can speed up load times
#'   (although not enormously, as the entire file must still be read) and reduce
#'   the final object's size.
#' @param TIC Boolean. If TRUE, the TIC is extracted rather than the BPC.
#'
#' @return A `data.table` with columns for retention time (rt), and intensity
#'   (int).
grabMzmlBPC <- function(xml_data, rtrange, TIC=FALSE, incl_polarity){
  ms1_xpath <- paste0('//d1:spectrum[d1:cvParam[@name="ms level" and ',
                      '@value="1"]][d1:cvParam[@name="base peak intensity"]]')

  ms1_nodes <- xml2::xml_find_all(xml_data, ms1_xpath)

  rt_vals <- grabSpectraRt(ms1_nodes)
  if(incl_polarity){pol_vals <- grabSpectraPolarity(ms1_nodes)}
  if(!is.null(rtrange)){
    ms1_nodes <- ms1_nodes[rt_vals%between%rtrange]
    rt_vals <- rt_vals[rt_vals%between%rtrange]
    if(incl_polarity){
      pol_vals <- pol_vals[pol_vals%between%rtrange]
    }
  }

  int_xpath <- ifelse(TIC, "total ion current", "base peak intensity")
  int_xpath_full <- paste0('d1:cvParam[@name="', int_xpath, '"]')
  int_nodes <- xml2::xml_find_all(ms1_nodes, xpath = int_xpath_full)
  int_vals <- as.numeric(xml2::xml_attr(int_nodes, "value"))
  if(incl_polarity){
    return(data.table(rt=rt_vals, int=int_vals, polarity=pol_vals))
  }
  return(data.table(rt=rt_vals, int=int_vals))
}

#' Extract the DAD data from an mzML nodeset
#'
#' @param xml_data An `xml2` nodeset, usually created by applying `read_xml` to
#'   an mzML file.
#' @param rtrange A vector of length 2 containing an upper and lower bound on
#'   retention times of interest. Providing a range here can speed up load times
#'   (although not enormously, as the entire file must still be read) and reduce
#'   the final object's size.
#' @param file_metadata Information about the file used to decode the binary
#'   arrays containing m/z and intensity information.
#' @return A `data.table` with columns for retention time (rt), wavelength
#' (lambda), and intensity (int).
grabMzmlDAD <- function(xml_data, rtrange, file_metadata){
  dad_xpath <- "//d1:spectrum[contains(@id,'controllerType=4')]"
  dad_nodes <- xml2::xml_find_all(xml_data, dad_xpath)
  if(!length(dad_nodes)){
    return(data.table(rt=numeric(), lambda=numeric(), int=numeric()))
  }

  rt_vals <- grabSpectraRt(dad_nodes)
  if(!is.null(rtrange)){
    dad_nodes <- dad_nodes[rt_vals%between%rtrange]
    rt_vals <- rt_vals[rt_vals%between%rtrange]
  }

  uv_vals <- grabSpectraMz(dad_nodes, file_metadata)
  int_vals <- grabSpectraInt(dad_nodes, file_metadata)

  int <- NULL #To prevent R CMD check "notes"  when using data.table syntax
  all_data <- data.table(rt=rep(rt_vals, sapply(uv_vals, length)),
                         lambda=unlist(uv_vals), int=as.numeric(unlist(int_vals)))
}


# Get spectrum things (functions of xml_nodes) ----

#' Extract the retention time from the spectra of an mzML nodeset
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#'   by the mass spectrometer, usually produced by applying `xml_find_all` to an
#'   MS1 or MS2 nodeset.
#'
#' @return A numeric vector of retention times, one for each scan
grabSpectraRt <- function(xml_nodes){
  rt_xpath <- 'd1:scanList/d1:scan/d1:cvParam[@name="scan start time"]'
  rt_nodes <- xml2::xml_find_all(xml_nodes, rt_xpath)
  rt_unit <- unique(xml_attr(rt_nodes, "unitName"))
  rt_vals <- as.numeric(xml2::xml_attr(rt_nodes, "value"))

  if(!"minute"%in%rt_unit) rt_vals=rt_vals/60
  # if(any(rt_vals>150)){
  #   # Guess RT is in seconds if the run is more than 150 long
  #   # A 2.5 minute run is unheard of, and a 2.5 hour run is unheard of
  #   rt_vals <- rt_vals/60
  # }

  rt_vals
}


#' Extract the precursor mass from the spectra of an mzML nodeset
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#'   by the mass spectrometer, usually produced by applying `xml_find_all` to an
#'   MS1 or MS2 nodeset.
#'
#' @return A numeric vector of precursor masses, one for each scan
grabSpectraPremz <- function(xml_nodes){
  premz_xpath <- paste0('d1:precursorList/d1:precursor/d1:selectedIonList',
                        '/d1:selectedIon/d1:cvParam[@name="selected ion m/z"]')
  premz_nodes <- xml2::xml_find_all(xml_nodes, premz_xpath)
  as.numeric(xml2::xml_attr(premz_nodes, "value"))
}


#' Extract the collison energies from the spectra of an mzML nodeset
#'
#' Although the collision energy is typically fixed per file, it's equally fast
#' (afaik) to just grab them all individually here. Also, I'm worried about
#' these rumors of "ramped" collision energies
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#'   by the mass spectrometer, usually produced by applying `xml_find_all` to an
#'   MS1 or MS2 nodeset.
#'
#' @return A numeric vector of collision energies, one for each scan.
grabSpectraVoltage <- function(xml_nodes){
  volt_xpath <- paste0('d1:precursorList/d1:precursor/d1:activation',
                       '/d1:cvParam[@name="collision energy"]')
  volt_nodes <- xml2::xml_find_all(xml_nodes, volt_xpath)
  if(length(volt_nodes)==0)return(NA_integer_)
  as.integer(xml2::xml_attr(volt_nodes, "value"))
}

grabSpectraPolarity <- function(xml_nodes){
  pol_xpath <- 'd1:cvParam[@name="positive scan" or @name="negative scan"]'
  pol_nodes <- xml2::xml_find_all(xml_nodes, pol_xpath)
  pol_vals <- xml2::xml_attr(pol_nodes, "name")
  num_pol <- numeric(length(pol_nodes))
  num_pol[pol_vals=="positive scan"] <- 1
  num_pol[pol_vals=="negative scan"] <- -1
  num_pol
}



#' Extract the mass-to-charge data from the spectra of an mzML nodeset
#'
#' The mz and intensity information of mzML files are encoded as binary arrays,
#' sometimes compressed via gzip or zlib or numpress. This code finds all the
#' m/z binary arrays and converts them back to the original measurements. See
#' https://github.com/ProteoWizard/pwiz/issues/1301
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#'   by the mass spectrometer, usually produced by applying `xml_find_all` to an
#'   MS1 or MS2 nodeset.
#' @param file_metadata Information about the file used to decode the binary
#'   arrays containing m/z and intensity information. Here, the compression and
#'   mz precision information is relevant.
#'
#' @return A numeric vector of masses, many for each scan.
grabSpectraMz <- function(xml_nodes, file_metadata){
  mz_xpath <- 'd1:binaryDataArrayList/d1:binaryDataArray[1]/d1:binary'
  mz_vals <- xml2::xml_text(xml2::xml_find_all(xml_nodes, mz_xpath))
  lapply(mz_vals, function(binary){
    if(!nchar(binary))return(numeric(0))
    decoded_binary <- base64enc::base64decode(binary)
    raw_binary <- as.raw(decoded_binary)
    decomp_binary <- memDecompress(raw_binary, type = file_metadata$compression)
    final_binary <- readBin(decomp_binary, what = "double",
                            n=length(decomp_binary)/file_metadata$mz_precision,
                            size = file_metadata$mz_precision)
  })
}


#' Extract the intensity information from the spectra of an mzML nodeset
#'
#' The mz and intensity information of mzML files are encoded as binary arrays,
#' sometimes compressed via gzip or zlib or numpress. This code finds all the
#' intensity binary arrays and converts them back to the original measurements.
#' See https://github.com/ProteoWizard/pwiz/issues/1301
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#'   by the mass spectrometer, usually produced by applying `xml_find_all` to an
#'   MS1 or MS2 nodeset.
#' @param file_metadata Information about the file used to decode the binary
#'   arrays containing m/z and intensity information. Here, the compression and
#'   int precision information is relevant.
#'
#' @return A numeric vector of intensities, many for each scan.
grabSpectraInt <- function(xml_nodes, file_metadata){
  int_xpath <- 'd1:binaryDataArrayList/d1:binaryDataArray[2]/d1:binary'
  int_vals <- xml2::xml_text(xml2::xml_find_all(xml_nodes, int_xpath))
  int_vals <- lapply(int_vals, function(binary){
    if(!nchar(binary))return(numeric(0))
    decoded_binary <- base64enc::base64decode(binary)
    raw_binary <- as.raw(decoded_binary)
    decomp_binary <- memDecompress(raw_binary, type = file_metadata$compression)
    final_binary <- readBin(decomp_binary, what = "double",
                            n=length(decomp_binary)/file_metadata$int_precision,
                            size = file_metadata$int_precision)
  })
}



# Get chromatogram things (functions of xml_nodes) ----
grabMzmlChroms <- function(xml_data, file_metadata){
  chrom_xpath <- '//d1:chromatogram'
  chrom_nodes <- xml2::xml_find_all(xml_data, chrom_xpath)

  chrom_id <- xml_attr(chrom_nodes, "id")
  chrom_idx <- xml_attr(chrom_nodes, "index")
  target_mz_xpath <- 'd1:precursor//d1:cvParam[@name="isolation window target m/z"]'
  target_mzs <- as.numeric(xml_attr(xml_child(chrom_nodes, target_mz_xpath), "value"))
  product_mz_xpath <- 'd1:product//d1:cvParam[@name="isolation window target m/z"]'
  product_mzs <- as.numeric(xml_attr(xml_child(chrom_nodes, product_mz_xpath), "value"))

  time_vals <- grabChromRt(chrom_nodes, file_metadata)
  int_vals <- grabChromInt(chrom_nodes, file_metadata)

  all_data <- data.table(chrom_type=rep(chrom_id, lengths(time_vals)),
                         chrom_index=rep(chrom_idx, lengths(time_vals)),
                         target_mz=rep(target_mzs, lengths(time_vals)),
                         product_mz=rep(product_mzs, lengths(time_vals)),
                         rt=unlist(time_vals), int=as.numeric(unlist(int_vals)))
}

grabChromRt <- function(chrom_nodes, file_metadata){
  # Exactly the same as grabbing the mz values for spectra
  # In chromatograms, the first vector is time
  # In spectra, the first vector is mz
  # We index by the first vector so these are identical
  grabSpectraMz(chrom_nodes, file_metadata)
}
grabChromInt <- function(chrom_nodes, file_metadata){
  # Exactly the same as for spectra
  # Second chromatogram is intensity in both cases
  grabSpectraInt(chrom_nodes, file_metadata)
}



# Other helper functions ----
shrinkRTrangemzML <- function(xml_nodes, rtrange){
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
