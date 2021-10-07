# TO-DO:
# Write documentation
# Write tests
# Write mzXML functions
# Publish new version to CRAN


# minifyMzml ----

#' Shrink mzML files by including only data points near masses of interest
#'
#' mzML files can be annoyingly large if only a few masses are of interest. This large size makes it
#' difficult to share them online for debugging purposes and often means that untargeted algorithms
#' spend a lot of time picking peaks in data that's irrelevant. minifyMzml is a function designed to
#' "minify" mzML files by extracting only those data points that are within a ppm error of an m/z value
#' of interest, and returns the file essentially otherwise unchanged. This function currently works
#' only on MS1 data, but is reasonably expandable if demand becomes evident.
#'
#' @param filename The name of a single file to be minified, usually produced by Proteowizard's `msconvert`
#' or something similar.
#' @param output_filename The name of the file to be written out.
#' @param mz_blacklist A vector of m/z values that should be excluded from the minified file. This argument
#' must be used with the `ppm` argument and should not be used with mz_whitelist. For each mass provided, an
#' m/z window of +/- `ppm` is calculated, and all data points within that window are removed.
#' @param mz_whitelist A vector of m/z values that should be included in the minified file. This argument
#' must be used with the `ppm` argument and should not be used with mz_blacklist. For each mass provided, an
#' m/z window of +/- `ppm` is calculated, and all data points within that window are kept.
#' @param ppm The parts-per-million error of the instrument used to collect the original file.
#' @param warn Boolean. Should the function warn the user when removing an index from an mzML file?
#'
#' @return Invisibly, the name of the new file.
#' @export
#'
#' @examples
#' \dontrun{
#' library(RaMS)
#' # Extract data corresponding to only valine and homarine
#' # m/z = 118.0865 and 138.0555, respectively
#' filename <- system.file("extdata", "LB12HL_AB.mzML.gz", package = "RaMS")
#' output_filename <- "mini_LB12HL_AB.mzML"
#' include_mzs <- c(118.0865, 138.0555)
#' minifyMzml(filename, output_filename, mz_whitelist=include_mzs, ppm=5)
#' unlink(output_filename)
#'
#' # Exclude data corresponding to valine and homarine
#' filename <- system.file("extdata", "LB12HL_AB.mzML.gz", package = "RaMS")
#' output_filename <- "mini_LB12HL_AB.mzML"
#' exclude_mzs <- c(118.0865, 138.0555)
#' minifyMzml(filename, output_filename, mz_blacklist=exclude_mzs, ppm=5)
#' unlink(output_filename)
#' }
minifyMzml <- function(filename, output_filename,
                         mz_blacklist=NULL, mz_whitelist=NULL,
                         ppm=NULL, warn=TRUE){
  xml_data <- xml2::read_xml(filename)

  RaMS:::checkFileType(xml_data, "mzML")
  file_metadata <- RaMS:::grabMzmlEncodingData(xml_data)

  # Check for indexed mzML and drop index if present, with warning
  if(xml2::xml_name(xml_data)=="indexedmzML" && warn){
    warning(paste0("mzML file ", basename(filename), " contains an index. ",
                   "I don't know how to recompile indices so it's ",
                   "getting dropped for the minified file."))
    xml_data <- xml2::xml_new_root(xml2::xml_find_first(xml_data, "//d1:mzML"))
  }

  # Find MS1 intensity and m/z nodes
  ms1_xpath <- paste0('//d1:spectrum[d1:cvParam[@name="ms level" and ',
                      '@value="1"]][d1:cvParam[@name="base peak intensity"]]')
  ms1_nodes <- xml2::xml_find_all(xml_data, ms1_xpath)
  mz_xpath <- 'd1:binaryDataArrayList/d1:binaryDataArray[1]/d1:binary//text()'
  ms1_mz_nodes <- xml2::xml_find_all(ms1_nodes, mz_xpath)
  int_xpath <- 'd1:binaryDataArrayList/d1:binaryDataArray[2]/d1:binary//text()'
  ms1_int_nodes <- xml2::xml_find_all(ms1_nodes, int_xpath)

  # Convert MS1 nodes into data.tables
  ms1_minified <- mapply(function(ms1_mz_node, ms1_int_node){
    mzs <- getEncoded(xml2::xml_text(ms1_mz_node),
                      compression_type = file_metadata$compression,
                      bin_precision = file_metadata$mz_precision)
    ints <- getEncoded(xml2::xml_text(ms1_int_node),
                       compression_type = file_metadata$compression,
                       bin_precision = file_metadata$int_precision)
    if(!is.null(mz_whitelist)){
      whitelist_data <- lapply(mz_whitelist, function(mz_i){
        range_i <- pmppm(mz_i, ppm)
        mz_idxs <- mzs>min(range_i)&mzs<max(range_i)
        ints <- ints[mz_idxs]
        mzs <- mzs[mz_idxs]
        cbind(mzs=mzs, ints=ints)
      })
      if(all(sapply(whitelist_data, length)==0)){
        output_mat <- cbind(mzs=numeric(0), ints=numeric(0))
      } else {
        output_mat <- do.call(what = "rbind", whitelist_data)
      }
    } else if(!is.null(mz_blacklist)){
      iterated_mzs <- mzs
      iterated_ints <- ints
      for(mz_i in unique(mz_blacklist)){
        range_i <- pmppm(mz_i, ppm)
        iterate_idxs <- iterated_mzs>min(range_i)&iterated_mzs<max(range_i)
        iterated_mzs <- iterated_mzs[!iterate_idxs]
        iterated_ints <- iterated_ints[!iterate_idxs]
      }
      output_mat <- cbind(mzs=iterated_mzs, ints=iterated_ints)
      if(any(is.na(output_mat))){
        browser("Found you an NA:")
      }
    } else {
      stop("Either `mz_whitelist` or `mz_blacklist` must not be NULL")
    }
    if(nrow(output_mat)==0){
      recoded_mzs <- ""
      recoded_ints <- ""
      bpmz <- 0
      bpint <- 0
      ticur <- 0
      minmz <- 0
      maxmz <- 0

      arraylength <- 0
    } else {
      recoded_mzs <- giveEncoding(output_mat[,"mzs"],
                                  compression_type = file_metadata$compression,
                                  bin_precision = file_metadata$mz_precision)
      recoded_ints <- giveEncoding(output_mat[,"ints"],
                                   compression_type = file_metadata$compression,
                                   bin_precision = file_metadata$int_precision)
      bpmz <- unname(output_mat[which.max(output_mat[,"ints"]), "mzs"])
      bpint <- max(output_mat[,"ints"], na.rm = TRUE)
      ticur <- sum(output_mat[,"ints"], na.rm = TRUE)
      minmz <- min(output_mat[,"mzs"], na.rm = TRUE)
      maxmz <- max(output_mat[,"mzs"], na.rm = TRUE)

      arraylength <- nrow(output_mat)
    }
    c(mzs=recoded_mzs, ints=recoded_ints, bpmz=bpmz, bpint=bpint,
      ticur=ticur, minmz=minmz, maxmz=maxmz,
      mz_enclength=nchar(recoded_mzs),
      int_enclength=nchar(recoded_ints),
      arraylength=arraylength)
  }, ms1_mz_nodes, ms1_int_nodes)
  ms1_minified <- as.data.frame(t(ms1_minified))

  xml2::xml_text(ms1_mz_nodes) <- ms1_minified$mzs
  xml2::xml_text(ms1_int_nodes) <- ms1_minified$ints

  xml2::xml_attr(ms1_nodes, "defaultArrayLength") <- ms1_minified$arraylength

  mz_enclength_xpath <- "d1:binaryDataArrayList/d1:binaryDataArray[d1:cvParam[@name='m/z array']]"
  mz_enclength_nodes <- xml2::xml_find_all(ms1_nodes, mz_enclength_xpath)
  xml2::xml_attr(mz_enclength_nodes, "encodedLength") <- ms1_minified$mz_enclength

  int_enclength_xpath <- "d1:binaryDataArrayList/d1:binaryDataArray[d1:cvParam[@name='intensity array']]"
  int_enclength_nodes <- xml2::xml_find_all(ms1_nodes, int_enclength_xpath)
  xml2::xml_attr(int_enclength_nodes, "encodedLength") <- ms1_minified$int_enclength


  bpmz_xpath <- "d1:cvParam[@name='base peak m/z']"
  init_bpmz_node <- xml2::xml_find_all(ms1_nodes, bpmz_xpath)
  xml2::xml_attr(init_bpmz_node, "value") <- ms1_minified$bpmz

  bpint_xpath <- "d1:cvParam[@name='base peak intensity']"
  init_bpc_node <- xml2::xml_find_all(ms1_nodes, bpint_xpath)
  xml2::xml_attr(init_bpc_node, "value") <- ms1_minified$bpint

  ticur_xpath <- "d1:cvParam[@name='total ion current']"
  init_tic_node <- xml2::xml_find_all(ms1_nodes, ticur_xpath)
  xml2::xml_attr(init_tic_node, "value") <- ms1_minified$ticur

  minmz_xpath <- "d1:cvParam[@name='lowest observed m/z']"
  init_minmz_node <- xml2::xml_find_all(ms1_nodes, minmz_xpath)
  xml2::xml_attr(init_minmz_node, "value") <- ms1_minified$minmz

  maxmz_xpath <- "d1:cvParam[@name='highest observed m/z']"
  init_maxmz_node <- xml2::xml_find_all(ms1_nodes, maxmz_xpath)
  xml2::xml_attr(init_maxmz_node, "value") <- ms1_minified$maxmz


  # Add note that RaMS was used to shrink the file
  proclist_node <- xml2::xml_find_all(xml_data, "//d1:dataProcessingList")
  xml2::xml_add_child(proclist_node, "dataProcessing", id="RaMS_R_package")
  proc_node <- xml2::xml_find_all(proclist_node, '//dataProcessing[@id="RaMS_R_package"]')
  xml2::xml_add_child(proc_node, "processingMethod", order=0, softwareRef="RaMS")
  meth_node <- xml2::xml_find_all(proclist_node, '//processingMethod[@order="0"]')
  if(!is.null(mz_whitelist)){
    xml2::xml_add_child(meth_node, "userParam", cvRef="MS", accession="MS:1009000",
                        name="Minification by m/z whitelist",
                        value=paste0(mz_whitelist, collapse = "; "))
  } else {
    xml2::xml_add_child(meth_node, "userParam", cvRef="MS", accession="MS:1009001",
                        name="Minification by m/z blacklist",
                        value=paste0(mz_blacklist, collapse = "; "))
  }
  xml2::xml_add_child(meth_node, "userParam", cvRef="MS", accession="MS:1009002",
                      name="ppm error for minification", value=ppm)
  process_count <- as.numeric(xml2::xml_attr(proclist_node, "count"))+1
  xml2::xml_attr(proclist_node, "count") <- process_count
  softlist_node <- xml2::xml_find_all(xml_data, "//d1:softwareList")
  xml2::xml_add_child(softlist_node, "software", id="RaMS",
                      version=as.character(packageVersion("RaMS")))
  soft_node <- xml2::xml_find_all(softlist_node, 'software[@id="RaMS"]')
  xml2::xml_add_child(soft_node, "userParam", name="RaMS R package")
  software_count <- as.numeric(xml2::xml_attr(softlist_node, "count"))+1
  xml2::xml_attr(softlist_node, "count") <- software_count


  # And write out the new version
  xml2::write_xml(x = xml_data, file = output_filename)
  return(invisible(output_filename))
}




#' Convert from compressed binary to R numeric vector
#'
#' @param mzint_nodes The XML nodes containing the compressed binary string
#' @param compression_type Compression type to be used by memDecompress
#' @param bin_precision The bit (?) precision used by readBin
#'
#' @return A numeric vector of m/z or intensity values
getEncoded <- function(mzint_nodes, compression_type, bin_precision){
  decoded_mzs <- base64enc::base64decode(mzint_nodes)
  decomp_mzs <- memDecompress(decoded_mzs, type = compression_type)
  readBin(decomp_mzs, what = "double", n=length(decomp_mzs)/bin_precision,
          size = bin_precision)
}

#' Convert from R numeric vector to compressed binary
#'
#' @param mzint_vals A numeric vector of m/z or intensity values
#' @param compression_type Compression type to be used by memCompress
#' @param bin_precision The bit (?) precision used by writeBin
#'
#' @return A single base64-encoded string of compressed binary values
giveEncoding <- function(mzint_vals, compression_type, bin_precision){
  comp_ints <- writeBin(mzint_vals, raw(0), size = bin_precision)
  new_raw_ints <- memCompress(comp_ints, type=compression_type)
  base64enc::base64encode(new_raw_ints)
}






# Test it ----
library(tidyverse)
library(data.table)
library(RaMS)
ppm <- 5
filename <- list.files(r"(G:\My Drive\AllMeso\mzMLs\pos)", full.names = TRUE)[3]
output_filename <- "C:/Users/willi/Desktop/mini_mzML.mzML"
msdata_init <- RaMS::grabMSdata(filename)

### Whitelist ----
mz_whitelist <- c(118.0865, 138.0555)
minifyMzml(filename, output_filename = "C:/Users/willi/Desktop/mini_mzML.mzML",
             mz_blacklist=NULL, mz_whitelist=mz_whitelist, ppm=ppm)










msdata <- RaMS::grabMSdata(output_filename)

msdata$MS1[mz%between%pmppm(118.0865)] %>%
  ggplot() + geom_line(aes(x=rt, y=int)) + xlim(7, 8)
msdata_init$MS1[mz%between%pmppm(118.0865)] %>%
  ggplot() + geom_line(aes(x=rt, y=int)) + xlim(7, 8)

msdata$MS1[mz%between%pmppm(138.0555, 5)] %>%
  ggplot() + geom_line(aes(x=rt, y=int))
msdata_init$MS1[mz%between%pmppm(138.0555, 5)] %>%
  ggplot() + geom_line(aes(x=rt, y=int))
msdata$BPC %>%
  ggplot() + geom_line(aes(x=rt, y=int))
msdata_init$BPC %>%
  ggplot() + geom_line(aes(x=rt, y=int))
# Make sure all the data survived, only thing changing is filename
all.equal(msdata_init$MS1[mz%between%pmppm(138.0555)],
          msdata$MS1[mz%between%pmppm(138.0555)])

# Nothing there!
msdata$MS1[mz%between%pmppm(179.005326)]
# Unlike before!
msdata_init$MS1[mz%between%pmppm(179.005326)]


### Blacklist ----
mz_blacklist <- c(138.0555, 118.0865)
minifyMzml(filename, output_filename = "C:/Users/willi/Desktop/mini_mzML.mzML",
             mz_blacklist=mz_blacklist, mz_whitelist=NULL, ppm=ppm)


msdata <- grabMSdata(output_filename)

msdata$MS1[mz%between%pmppm(118.0865, 5)] %>%
  ggplot() + geom_line(aes(x=rt, y=int))
msdata_init$MS1[mz%between%pmppm(118.0865, 5)] %>%
  ggplot() + geom_line(aes(x=rt, y=int))

msdata$MS1[mz%between%pmppm(138.0555, 5)] %>%
  ggplot() + geom_line(aes(x=rt, y=int))
msdata_init$MS1[mz%between%pmppm(138.0555, 5)] %>%
  ggplot() + geom_line(aes(x=rt, y=int))
msdata$BPC %>%
  ggplot() + geom_line(aes(x=rt, y=int))
msdata_init$BPC %>%
  ggplot() + geom_line(aes(x=rt, y=int))

all.equal(msdata_init$MS1[mz%between%pmppm(138.0555)],
          msdata$MS1[mz%between%pmppm(138.0555)])




### XCMS check ----
library(xcms)
raw_data <- readMSData(output_filename, pdata = NULL, msLevel. = 1)
bpis <- chromatogram(raw_data, aggregationFun = "max")
plot(bpis)
