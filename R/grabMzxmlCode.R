
#' @import xml2
#' @importFrom base64enc base64decode
NULL

#' Grab the BPC or TIC for a given mzXML file
#'
#' @param filename The name of the mzXML file to be read.
#' @param TIC Should the BPC or TIC be read? If TIC=TRUE,
#' the total ion current (TIC)
#' for each retention time is read in. Otherwise, the BPC
#' (base peak chromatogram) is read for each retention time, corresponding
#' to the *maximum* intensity for each scan.
#'
#' @return A data.frame object with columns for retention time (rt) in minutes and
#' and intensity (int) corresponding to the TIC or BPC for a given file.
#'
#' @export
grabMzxmlBPC <- function(filename, TIC=FALSE){
  mz_xml <- xml2::read_xml(filename)

  scan_nodes <- xml2::xml_find_all(mz_xml, '//d1:scan')
  rt_chrs <- xml2::xml_attr(scan_nodes, "retentionTime")
  rt_vals <- as.numeric(gsub(pattern = "PT|S", replacement = "", rt_chrs))/60

  int_attr <- ifelse(TIC, "totIonCurrent", "basePeakIntensity")
  int_vals <- as.numeric(xml2::xml_attr(scan_nodes, int_attr))

  return(data.frame(rt=rt_vals, int=int_vals))
}



#' Read an mzXML file into a data.frame
#'
#' @details This function reads an mzXML file into R's working memory. mzXML files
#' are fundamentally XML documents, which allows rapid access to the data by
#' parsing the XML. The R package `xml2::` is used for this purpose here.
#' Retention time information can be read directly, while *m/z* and intensity
#' information must be decoded from binary. To read data from the mzML format
#' instead, check out grabMzmlData.
#'
#' @param filename The name of the mzXML file to be read.
#'
#' @return A data.frame object with columns for retention time (rt) in minutes,
#' m/z (mz), and intensity (int).
#'
#' @export
grabMzxmlData <- function(filename){
  xml_data <- xml2::read_xml(filename)

  rt_xpath <- '//d1:scan'
  rt_nodes <- xml2::xml_find_all(xml_data, rt_xpath)
  rt_attrs <- xml2::xml_attr(rt_nodes, "retentionTime")
  rt_vals <- as.numeric(gsub("PT|S", "", rt_attrs))

  peak_metadata <- grabMzxmlMetadata(xml_data)

  all_peak_nodes <- xml2::xml_text(xml2::xml_find_all(xml_data, '//d1:peaks'))
  vals <- lapply(all_peak_nodes, function(binary){
    decoded_binary <- base64enc::base64decode(binary)
    raw_binary <- as.raw(decoded_binary)
    decomp_binary <- memDecompress(raw_binary, type = peak_metadata$compression)
    final_binary <- readBin(decomp_binary, what = "numeric",
                            n=length(decomp_binary)/peak_metadata$precision,
                            size = peak_metadata$precision,
                            endian = peak_metadata$endi_enc)
    matrix(final_binary, ncol = 2, byrow = TRUE)
  })
  add_rts <- mapply(cbind, rt_vals, vals)
  output <- as.data.frame(do.call(add_rts, what = "rbind")%*%diag(c(1/60, 1, 1)))
  names(output) <- c("rt", "mz", "int")
  output
}



#' Helper function to extract mzXML file metadata
#'
#' @param xml_data mzXML data parsed by xml2
#'
#' @return A list of values used by other parsing functions, such as
#' compression type and precision
#'
#' @export
grabMzxmlMetadata <- function(xml_data){
  peak_metadata <- xml2::xml_find_first(xml_data, '//d1:peaks')
  compr_type <- xml2::xml_attr(peak_metadata, "compressionType")
  compr <- switch(compr_type,
                  `zlib compression`="gzip",
                  `none`="none")

  enc_type <- xml2::xml_attr(peak_metadata, "precision")
  precision <- as.numeric(enc_type)/8

  byte_order <- xml2::xml_attr(peak_metadata, "byteOrder")
  endi_enc <- switch(byte_order, `network`="big")

  list(compression=compr, precision=precision, endi_enc=endi_enc)
}

