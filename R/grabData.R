# Requires xml2, mzR, data.table

#' Read an mzML file into a data.frame
#'
#' @param filename The name of the mzML file to be read
#'
#' @return A data.frame object with columns for retention time (rt), m/z (mz),
#' and intensity (int)
#' @export
#'
#' @examples
#' mzML_filename <- "data/180205_Poo_TruePoo_Full2.mzML"
#' grabSingleFileData(mzML_filename)
grabSingleFileData <- function(filename){
  msdata <- mzR:::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  spectra_list <- lapply(seq_len(nrow(fullhd)), function(x){
    given_peaks <- mzR::peaks(msdata, x)
    rtime <- fullhd[x, "retentionTime"]
    return(cbind(rtime, given_peaks))
  })
  all_data <- as.data.table(do.call(rbind, spectra_list))
  names(all_data) <- c("rt", "mz", "int")
  return(all_data)
}



#' Grab the BPC or TIC for a given file
#'
#' @param filename The name of the mzML file to be read
#' @param TIC Should the BPC or TIC be read? If TIC=TRUE,
#' the total ion current (TIC)
#' for each retention time is read in. Otherwise, the BPC
#' (base peak chromatogram) is read for each retention time, corresponding
#' to the *maximum* intensity for each scan.
#'
#' @return A data.frame object with columns for retention time (rt) and
#' and intensity (int) corresponding to the
#' @export
#'
#' @examples
#' mzML_filename <- "data/180205_Poo_TruePoo_Full2.mzML"
#' grabSingleFileBPC(mzML_filename)
#' grabSingleFileBPC(mzML_filename, TIC=TRUE)
grabSingleFileBPC <- function(filename, TIC=FALSE){
  mz_xml <- xml2::read_xml(filename)

  rt_nodes <- xml2::xml_find_all(mz_xml, '//d1:cvParam[@name="scan start time"]')
  rt_vals <- as.numeric(xml2::xml_attr(rt_nodes, "value"))

  int_xpath <- ifelse(TIC, "total ion current", "base peak intensity")
  int_xpath_full <- paste0('//d1:cvParam[@name="', int_xpath, '"]')
  int_nodes <- xml2::xml_find_all(mz_xml, xpath = int_xpath_full)
  int_vals <- xml2::xml_attr(int_nodes, "value")
  return(data.table)
}

