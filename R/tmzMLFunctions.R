# TO-DO:
# Add glimpse() functionality
# Include TIC/BPC as separate node
# Fix MS2 voltage encoding issues
# Finish documentation
# Add warnings/error handling
# Add tests
# Add .tmzML minified file for demos and examples

addEncNode <- function(parent_node, dubset, name){
  new_node <- xml2::xml_add_child(parent_node, name)
  xml2::xml_text(new_node) <- giveEncoding(
    dubset[[name]], compression_type = "gzip",
    bin_precision = 8, endi_enc = "little"
  )
}
addMS1 <- function(dubset, ms1_node){
  node_i <- xml2::xml_add_child(ms1_node, "dubset")
  addEncNode(node_i, dubset, "rt")
  addEncNode(node_i, dubset, "mz")
  addEncNode(node_i, dubset, "int")
  xml2::xml_attr(node_i, "minmz") <- min(dubset$mz)
  xml2::xml_attr(node_i, "maxmz") <- max(dubset$mz)
  return(NULL)
}
addMS2 <- function(dubset, ms2_node){
  node_i <- xml2::xml_add_child(ms2_node, "dubset")
  addEncNode(node_i, dubset, "rt")
  addEncNode(node_i, dubset, "premz")
  addEncNode(node_i, dubset, "fragmz")
  addEncNode(node_i, dubset, "voltage")
  addEncNode(node_i, dubset, "int")
  xml2::xml_attr(node_i, "minmz") <- min(dubset$premz)
  xml2::xml_attr(node_i, "maxmz") <- max(dubset$premz)
  return(NULL)
}

#' Maker of tmzML documents
#'
#' This function converts mzML and mzXML documents into "transposed" mzML
#' (tmzML) documents. Traditional mass-spec data is organized by scan number,
#' corresponding to retention time, but this isn't always the most sensible
#' format. Often, it makes more sense to organize a mass-spec file by m/z ratio
#' instead. This allows parsers to scan and decode a much smaller portion of
#' the file when searching for a specific mass, as opposed to the traditional
#' format which requires that every scan be opened, searched, and subset. The
#' tmzML document implements this strategy and allows the creation of MS object
#' representations that use essentially zero memory because the data is read
#' off the disk instead of being stored in RAM. RaMS has been designed to
#' interface with these new file types identically to traditional files,
#' allowing all your favorite tidyverse tricks to work just as well and much
#' more quickly.
#'
#' @param input_filename Character vector of length 1 with the name of the file
#' to be converted. Can only handle
#' mzML and mzXML currently - other formats should be converted to one of
#' these first, using (for example) Proteowizard's msconvert tool.
#' @param output_filename The name of the file that will be written out. Should
#' end in ".tmzML" and will throw a warning otherwise. Often, it makes sense
#' to have two folders in a working directory, one containing the original
#' mzML files and a second, parallel folder for the tmzMLs.
#' @param verbosity Numeric value between 0 and 2, corresponding to level of
#' verbosity shared by the function as it proceeds. 0 means no output, 1 will
#' produce mile markers after file opening, MS1 and MS2 conversion, and 2 will
#' provide progress bars between each mile marker.
#' @param binwidth Numeric value controlling the width of the bins in m/z space to create. Because MS data
#' is created in such a way that m/z values are continuous, they must be binned
#' together to create a discrete representation that can be searched efficiently.
#' Lower values (0.1-1) will have faster retrieval times, while higher values (5-10)
#' will have faster conversion times.
#'
#' @return An msdata_connection object. This object behaves exactly like a normal
#' RaMS list with values for MS1, MS2, etc. but secretly just contains pointers
#' to the files requested because the data is extracted on the fly. The S3
#' msdata_connection object is necessary to create new behaviors for `$` and `[` that
#' allow indexing like normal.
#' @export
#'
#' @examples
#' \dontrun{
#' sample_dir <- system.file("extdata", package = "RaMS")
#' sample_files <- list.files(sample_dir, full.names=TRUE, pattern="LB.*mzML")
#' tmzml_filenames <- gsub(x=sample_files, "\\.mzML.gz", ".tmzML")
#'
#' # Convert a single file
#' tmzmlMaker(sample_files[1], tmzml_filenames[1])
#' file_data <- grabMSdata(tmzml_filenames[1], grab_what="everything", verbosity=2)
#' file_data$MS1[mz%between%pmppm(118.0865)]
#'
#' # Multiple files
#' mapply(tmzmlMaker, sample_files, tmzml_filenames)
#' file_data <- grabMSdata(tmzml_filenames, grab_what="everything", verbosity=2)
#' betaine_data <- file_data$MS1[mz%between%pmppm(118.0865)]
#'
#' # Plot output
#' plot(betaine_data$rt, betaine_data$int, type="l")
#' library(ggplot2)
#' ggplot(betaine_data) + geom_line(aes(x=rt, y=int, color=filename))
#'
#' # Clean up afterward
#' file.remove(tmzml_filenames)
#' }
tmzmlMaker <- function(input_filename, output_filename=NULL,
                       verbosity=0, binwidth=3){
  if(is.null(output_filename)){
    output_filename <- gsub("\\.mzX?ML", ".tmzML", input_filename)
  }
  if(length(input_filename)>1){
    stop(paste("tmzmlMaker is not vectorized and cannot handle multiple",
               "files at once. Consider wrapping in 'mapply' instead."))
  }
  if(!endsWith(output_filename, "tmzML")){
    warning("The provided output_filename does not end in tmzML.")
  }

  if(verbosity>0){
    message("Reading in original file...")
  }
  msdata <- grabMSdata(input_filename, verbosity = 0)

  tmz_doc <- xml2::xml_new_root("tmzML")
  glimpse_node <- xml2::xml_add_child(tmz_doc, "glimpse")
  data_node <- xml2::xml_add_child(tmz_doc, "data")

  # Create glimpse ----
  xml2::xml_add_child(glimpse_node, "MS1")
  xml2::xml_add_child(glimpse_node, "MS2")

  # Handle MS1 data ----
  ms1_node <- xml2::xml_add_child(data_node, "MS1")
  msdata$MS1$cat <- floor(msdata$MS1$mz/binwidth)
  if(verbosity>0){
    message("Indexing MS1 data by m/z...")
  }
  split_MS1_mzs <- split(msdata$MS1, msdata$MS1$cat)
  if(verbosity>1){
    pb <- txtProgressBar(max = length(split_MS1_mzs), style = 3)
    for(i in seq_along(split_MS1_mzs)){
      addMS1(split_MS1_mzs[[i]], ms1_node)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  } else {
    for(i in seq_along(split_MS1_mzs)){
      addMS1(split_MS1_mzs[[i]], ms1_node)
    }
  }

  # Handle MS2 data ----
  ms2_node <- xml2::xml_add_child(data_node, "MS2")
  msdata$MS2$cat <- floor(msdata$MS2$premz/binwidth)
  if(verbosity>0){
    message("Indexing MS2 data by m/z...")
  }
  split_MS2_mzs <- split(msdata$MS2, msdata$MS2$cat)
  if(verbosity>1 & nrow(msdata$MS2)>0){
    pb <- txtProgressBar(max = length(split_MS2_mzs), style = 3)
    for(i in seq_along(split_MS2_mzs)){
      addMS2(split_MS2_mzs[[i]], ms2_node)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  } else {
    for(i in seq_along(split_MS2_mzs)){
      addMS2(split_MS2_mzs[[i]], ms2_node)
    }
  }

  # Write out ----
  xml2::write_xml(tmz_doc, output_filename)
  output_filename
}


#' Convert node to data.table
#'
#' @param dubset_node The "data subset" node with children rt, mz, etc.
#' @param ms_level The requested MS level to search for
#'
#' @return A data.table with columns depending on the MS level requested
node2dt <- function(dubset_node, ms_level){
  raw_encoded <- xml2::xml_text(xml2::xml_children(dubset_node))
  data_list <- lapply(raw_encoded, getEncoded, compression_type="gzip",
                      bin_precision = 8, endi_enc="little")
  if(ms_level=="MS1"){
    if(length(data_list)==0){
      data.table(rt=numeric(0), mz=numeric(0), int=numeric(0))
    } else if(length(data_list)==3){
      data.table(rt=data_list[[1]], mz=data_list[[2]], int=data_list[[3]])
    } else if(length(data_list)>3){
      rt_idxs <- seq(1, length(data_list), by=3)
      mz_idxs <- seq(2, length(data_list), by=3)
      int_idxs <- seq(3, length(data_list), by=3)
      disorg_table <- data.table(
        rt=unlist(data_list[rt_idxs]),
        mz=unlist(data_list[mz_idxs]),
        int=unlist(data_list[int_idxs])
      )
      mz <- NULL #To prevent R CMD check "notes"  when using data.table syntax
      rt <- NULL #To prevent R CMD check "notes"  when using data.table syntax

      disorg_table[order(rt, mz)]
    } else {
      stop("How???")
    }
  } else if(ms_level=="MS2"){
    if(length(data_list)==0){
      data.table(rt=numeric(0), premz=numeric(0), fragmz=numeric(0),
                 voltage=numeric(0), int=numeric(0))
    } else if(length(data_list)==5){
      data.table(rt=data_list[[1]], premz=data_list[[2]], fragmz=data_list[[3]],
                 voltage=data_list[[4]], int=data_list[[5]])
    } else if(length(data_list)>3){
      rt_idxs <- seq(1, length(data_list), by=5)
      premz_idxs <- seq(2, length(data_list), by=5)
      fragmz_idxs <- seq(3, length(data_list), by=5)
      voltage_idxs <- seq(4, length(data_list), by=5)
      int_idxs <- seq(5, length(data_list), by=5)
      disorg_table <- data.table(
        rt=unlist(data_list[rt_idxs]),
        premz=unlist(data_list[mz_idxs]),
        fragmz=unlist(data_list[mz_idxs]),
        voltage=unlist(data_list[mz_idxs]),
        int=unlist(data_list[int_idxs])
      )
      rt <- NULL #To prevent R CMD check "notes"  when using data.table syntax
      disorg_table[order(rt)]
    } else {
      stop("How???")
    }
  } else {
    stop(paste("MS level", ms_level, "not currently supported!"))
  }
}

#' S3 print option for msdata_connection objects
#'
#' @param x An msdata_connection object containing files and grab_what
#' @param ... Other arguments to be passed to print.default, I guess
#'
#' @return Messages, mostly
#' @export
"print.msdata_connection" <- function(x, ...){
  message("Hey, I'm not actually an object, sorry!")
  message("But you can pretend I'm a list containing data.tables:")
  message(paste(x[["connection"]][["grab_what"]], collapse = "; "))
  message("from the following files:")
  n_files <- length(x[["connection"]][["files"]])
  excess_file_message <- if(n_files>10){
    message(paste(head(x[["connection"]][["files"]]), collapse = "\n"))
    message(paste("and", n_files-6, "others"))
  } else {
    message(paste(x[["connection"]][["files"]], collapse = "\n"))
  }
  message("and access the data inside with $ and [ subsetting")
}

#' S3 dollar sign notation for msdata_connection objects
#'
#' @param msdata_obj An msdata_connection object containing files and grab_what
#' @param ms_level The requested MS level of the object
#'
#' @return An msdata_connection object with only a single MS level
#' @export
"$.msdata_connection" <- function(msdata_obj, ms_level){
  if(ms_level%in%msdata_obj[["connection"]][["grab_what"]]){
    msdata_obj[["connection"]][["grab_what"]] <- ms_level
  } else if(ms_level=="connection") {
    return(msdata_obj[["connection"]])
  } else {
    stop(paste0("It doesn't look like you requested '", ms_level,
                "' with grab_what when you created this object."))
  }
  return(msdata_obj)
}

#' S3 indexing for msdata_connection objects
#'
#' This is the step that actually performs the file opening and extraction!
#'
#' @param msdata_obj An msdata_connection object containing files and grab_what
#' @param sub_func The function that will be parsed and used to subset the file
#'
#' @return A data.table with columns rt, mz, int, and filename
#' @export
"[.msdata_connection" <- function(msdata_obj, sub_func){
  if(length(msdata_obj[["connection"]][["grab_what"]])!=1){
    stop("tmzML objects must first be indexed with $ to specify which data you'd like")
  }
  isub <- substitute(sub_func)
  function_name <- as.character(isub[[1]])
  col_name <- as.character(isub[[2]])
  if(!col_name%in%c("mz", "premz")){
    stop("tmzML documents currently only support subsetting by mz!")
  }
  mz_lims <- eval.parent(isub[[3]])
  if(length(mz_lims)==1){
    mz_lims <- c(eval.parent(isub[[3]]), eval.parent(isub[[4]]))
  }

  if(msdata_obj[["connection"]][["verbosity"]]>0){
    pb <- txtProgressBar(max = length(msdata_obj[["connection"]][["files"]]), style = 3)
  }
  allfile_list <- lapply(msdata_obj[["connection"]][["files"]], function(filename){
    tmzml <- xml2::read_xml(filename)

    mz_selectors <- paste0("@minmz<", max(mz_lims), " and @maxmz>", min(mz_lims))
    mz_xpath <- paste0("data/", msdata_obj[["connection"]][["grab_what"]],
                       "/dubset[", mz_selectors, "]")
    dubset_node <- xml2::xml_find_all(tmzml, mz_xpath)
    dubset_data <- node2dt(dubset_node, ms_level=msdata_obj[["connection"]][["grab_what"]])
    sub_data <- dubset_data[get(isub[[2]])%between%c(min(mz_lims), max(mz_lims))]
    sub_data$filename <- basename(filename)
    if(msdata_obj[["connection"]][["verbosity"]]>0){
      pb <- setTxtProgressBar(pb, which(msdata_obj[["connection"]][["files"]]==filename))
    }
    sub_data
  })
  if(msdata_obj[["connection"]][["verbosity"]]>0){
    close(pb)
  }
  rbindlist(allfile_list)
}


#' S3 constructor for msdata_connection
#'
#' @param x This is a thing?
#'
#' @return Itself, with the class?
#' @export
msdata_connection <- function(x){
  class(x) <- "msdata_connection"
  x
}
