# Requires xml2, base64enc

#' Grab the BPC or TIC for a given file
#'
#' @param filename The name of the mzML file to be read.
#' @param TIC Should the BPC or TIC be read? If TIC=TRUE,
#' the total ion current (TIC)
#' for each retention time is read in. Otherwise, the BPC
#' (base peak chromatogram) is read for each retention time, corresponding
#' to the *maximum* intensity for each scan.
#'
#' @return A data.frame object with columns for retention time (rt) in minutes and
#' and intensity (int) corresponding to the TIC or BPC for a given file.
#' @export
#'
#' @examples
#' mzML_filename <- system.file("extdata", "180205_Poo_TruePoo_Full2.mzML", package = "RaMS")
#' grabSingleFileBPC(mzML_filename)
#' grabSingleFileBPC(mzML_filename, TIC=TRUE)
grabSingleFileBPC <- function(filename, TIC=FALSE){
  mz_xml <- read_xml(filename)

  rt_nodes <- xml_find_all(mz_xml, '//d1:cvParam[@name="scan start time"]')
  rt_vals <- as.numeric(xml_attr(rt_nodes, "value"))

  int_xpath <- ifelse(TIC, "total ion current", "base peak intensity")
  int_xpath_full <- paste0('//d1:cvParam[@name="', int_xpath, '"]')
  int_nodes <- xml_find_all(mz_xml, xpath = int_xpath_full)
  int_vals <- xml_attr(int_nodes, "value")
  return(data.frame(rt=rt_vals, int=int_vals))
}



#' Read an mzML file into a data.frame
#'
#' @details This function reads an mzML file into R's working memory. mzML files
#' are fundamentally XML documents, which allows rapid access to the data by
#' parsing the XML. The R package `xml2::` is used for this purpose here.
#' Retention time information can be read directly, while *m/z* and intensity
#' information must be decoded from binary.
#'
#' @param filename The name of the mzML file to be read.
#'
#' @return A data.frame object with columns for retention time (rt) in minutes,
#' m/z (mz), and intensity (int).
#' @export
#'
#' @examples
#' mzML_filename <- system.file("extdata", "180205_Poo_TruePoo_Full2.mzML", package = "RaMS")
#' grabSingleFileData(mzML_filename)
grabSingleFileData <- function(filename){
  xml_data <- read_xml(filename)

  rt_nodes <- xml_find_all(xml_data, '//d1:cvParam[@name="scan start time"]')
  rt_vals <- as.numeric(xml_attr(rt_nodes, "value"))

  compression_xpath <- '//d1:cvParam[@accession="MS:1000574"]'
  compression_type <- xml_attr(xml_find_first(xml_data, compression_xpath), "name")
  compression <- switch(compression_type,
                        `zlib compression`="gzip",
                        `none`="none")

  mz_bit_xpath <- '//d1:cvParam[@accession="MS:1000523"]'
  mz_bit_type <- xml_attr(xml_find_first(xml_data, mz_bit_xpath), "name")
  mz_precision <- sub(mz_bit_type, pattern = "-bit float", replacement = "")
  mz_precision <- as.numeric(mz_precision)/8

  mz_xpath <- paste0('//d1:spectrum/d1:binaryDataArrayList',
                     '/d1:binaryDataArray[1]/d1:binary')
  mz_vals <- xml_text(xml_find_all(xml_data, mz_xpath))
  mz_vals <- lapply(mz_vals, function(binary){
      binary %>%
        base64enc::base64decode() %>%
        as.raw() %>%
        memDecompress(type = compression) %>%
        readBin(what = "double", n=length(.)/mz_precision,
                size = mz_precision)
    })

  int_bit_xpath <- '//d1:cvParam[@accession="MS:1000521"]'
  int_bit_type <- xml_attr(xml_find_first(xml_data, int_bit_xpath), "name")
  int_precision <- sub(int_bit_type, pattern = "-bit float", replacement = "")
  int_precision <- as.numeric(int_precision)/8

  int_xpath <- paste0('//d1:spectrum/d1:binaryDataArrayList',
                      '/d1:binaryDataArray[2]/d1:binary')
  int_vals <- xml_text(xml_find_all(xml_data, int_xpath))
  int_vals <- lapply(int_vals, function(binary){
    binary %>%
      base64enc::base64decode() %>%
      as.raw() %>%
      memDecompress(type = compression) %>%
      readBin(what = "double", n=length(.)/int_precision,
              size = int_precision)
  })

  data.frame(rt=rep(rt_vals, sapply(mz_vals, length)),
             mz=unlist(mz_vals), int=unlist(int_vals))
}
