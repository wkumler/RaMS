# TO-DO:
# Add glimpse() functionality
# Fix MS2 voltage encoding issues
# Finish documentation
# Add warnings/error handling
# Add tests
# Add tmzml minified file for demos and examples

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
tmzmlMaker <- function(input_filename, output_filename=NULL,
                       verbosity=0, binwidth=3){
  if(is.null(output_filename)){
    output_filename <- gsub("\\.mzX?ML", ".tmzML", input_filename)
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
  if(verbosity>1){
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
  message(paste(x[["grab_what"]], collapse = "; "))
  message("from the following files:")
  message(paste(x[["files"]], collapse = "\n"))
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
  msdata_obj$grab_what <- ms_level
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
  isub <- substitute(sub_func)
  function_name <- as.character(isub[[1]])
  col_name <- as.character(isub[[2]])
  if(!col_name%in%c("mz", "premz")){
    stop("tmzML documents currently only support subsetting by mz!")
  }
  mz_lims <- eval(isub[[3]])
  if(length(mz_lims)==1){
    mz_lims <- c(eval(isub[[3]]), eval(isub[[4]]))
  }

  allfile_list <- lapply(msdata_obj[["files"]], function(filename){
    tmzml <- xml2::read_xml(filename)

    mz_selectors <- paste0("@minmz<", max(mz_lims), " and @maxmz>", min(mz_lims))
    mz_xpath <- paste0("data/", msdata_obj[["grab_what"]],
                       "/dubset[", mz_selectors, "]")
    dubset_node <- xml2::xml_find_all(tmzml, mz_xpath)
    dubset_data <- node2dt(dubset_node, ms_level=msdata_obj[["grab_what"]])
    sub_data <- dubset_data[get(isub[[2]])%between%c(min(mz_lims), max(mz_lims))]
    sub_data$filename <- basename(filename)
    sub_data
  })
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
