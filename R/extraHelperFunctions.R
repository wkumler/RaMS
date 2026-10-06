
#' @export
data.table::`%between%`
#' @export
data.table::`%inrange%`

#' Plus/minus parts per million
#'
#' It shouldn't be hard to translate a point mass into a mass window bounded by
#' spectrometer accuracy.
#'
#' @param mass A length-1 numeric representing the mass of
#' interest for which a mass range is desired.
#' @param ppm The parts-per-million accuracy of the mass spectrometer on which
#' the data was collected.
#'
#' @return A length-2 numeric representing the mass range requested
#'
#' @seealso \code{\link{pmppm_vec}} for multiple masses at once
#'
#' @export
#'
#' @examples
#' pmppm(100, 5)
#' pmppm(1000000, 5)
#' pmppm(118.0865, 2.5)
#' pmppm(892.535313, 10)
pmppm <- function(mass, ppm=4){
  if(length(mass)>1)stop("Vectorized masses are not supported in pmppm; use pmppm_vec()")
  if(length(ppm)>1)stop("Vectorized ppms are not supported in pmppm")
  if(!is.numeric(mass))stop("'mass' must be numeric")
  if(!is.numeric(ppm))stop("'ppm' must be numeric")
  if(ppm<0)warning("'ppm' should be positive")
  c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))
}

#' Plus/minus parts per million for multiple masses
#'
#' A vectorized version of \code{\link{pmppm}} that converts any number of
#' masses into mass windows of the same ppm width. All the lower bounds are
#' returned together and all the upper bounds together, which is the format
#' expected by \code{\%inrange\%} and \code{\%between\%} (re-exported from
#' data.table) and by range joins.
#'
#' With \code{\%inrange\%}, each value is checked against every window, which
#' is useful for extracting data for several masses at once (e.g. inside
#' \code{dplyr::filter()}). With \code{\%between\%}, each value is checked only
#' against its own window, so \code{mass} should be the same length as the
#' values being checked.
#'
#' @section Use in tidyverse pipelines:
#' Because the bounds are grouped by type rather than by mass,
#' \code{pmppm_vec()} is not suited to creating a list-column with one window
#' per row: \code{mutate(mz_range=pmppm_vec(mz))} errors for most data
#' frames and, with exactly two rows, silently creates incorrect windows. For
#' one window per row, use \code{purrr::map(mz, pmppm)} instead, or add the
#' bounds as two columns with \code{mutate(as_tibble(pmppm_vec(mz)))}. See the
#' examples.
#'
#' @param mass A numeric vector of masses for which mass ranges are desired.
#' @param ppm A length-1 numeric giving the parts-per-million accuracy of the
#' mass spectrometer on which the data was collected. The same ppm is applied
#' to every mass.
#'
#' @return A list with two numeric elements, \code{lower} and \code{upper},
#' each the same length as \code{mass}, giving the bounds of the mass window
#' around each mass. A list is returned even when \code{mass} is length 1.
#'
#' @seealso \code{\link{pmppm}} for a single mass
#'
#' @export
#'
#' @examples
#' pmppm_vec(c(118.0865, 104.1070, 136.0618), 5)
#'
#' # Is each value within 5 ppm of any of the masses?
#' mz_vals <- c(104.1071, 104.2, 118.0866, 136.0619, 150)
#' mz_vals %inrange% pmppm_vec(c(118.0865, 104.1070, 136.0618), 5)
#'
#' \dontrun{
#' library(dplyr)
#' library(purrr)
#' library(ggplot2)
#' sample_dir <- system.file("extdata", package = "RaMS")
#' sample_files <- list.files(sample_dir, full.names=TRUE)
#' msdata <- grabMSdata(sample_files[c(3, 5, 6)], grab_what="MS1")
#'
#' cmpds <- tribble(
#'   ~cmpd_name, ~cmpd_mz,
#'   "Glycine betaine", 118.0865,
#'   "Choline", 104.1070,
#'   "Adenine", 136.0618
#' )
#'
#' # All data points within 5 ppm of any of the compounds
#' msdata$MS1 %>%
#'   filter(mz %inrange% pmppm_vec(cmpds$cmpd_mz, 5))
#'
#' # Label each data point with the compound it matched using a range join
#' cmpd_eics <- cmpds %>%
#'   mutate(as_tibble(pmppm_vec(cmpd_mz, 5))) %>%
#'   inner_join(msdata$MS1, join_by(between(y$mz, x$lower, x$upper)))
#' cmpd_eics %>%
#'   ggplot() +
#'   geom_line(aes(x=rt, y=int, color=filename)) +
#'   facet_wrap(~cmpd_name, ncol=1, scales="free_y")
#'
#' # For one window per row in a list-column, use pmppm() with map() rather
#' # than pmppm_vec(), which returns all the lower and upper bounds together
#' cmpds %>%
#'   mutate(mz_range=map(cmpd_mz, pmppm, ppm=5)) %>%
#'   mutate(eic=map(mz_range, ~filter(msdata$MS1, mz %between% .x)))
#'
#' # The same windows work with data.table syntax
#' msdata$MS1[mz %inrange% pmppm_vec(cmpds$cmpd_mz, 5)]
#' }
pmppm_vec <- function(mass, ppm=4){
  if(!is.numeric(mass))stop("'mass' must be numeric")
  if(!is.numeric(ppm) || length(ppm)!=1 || is.na(ppm))stop("'ppm' must be a length-1 numeric")
  if(ppm<0)warning("'ppm' should be positive")
  list(lower=mass*(1-ppm/1000000), upper=mass*(1+ppm/1000000))
}

#' Trapezoidal integration of mass-spec retention time / intensity values
#'
#' Performs a trapezoidal Riemann sum to calculate the area under the curve
#' for mass-spectrometry data. Accepts a vector of retention times and the
#' associated intensities and returns the area.
#'
#' @param rts A numeric vector of retention times across an MS peak. Should be
#' monotonically increasing and without duplicates or will throw a warning.
#' @param ints A numeric vector of measured intensities across an MS peak
#' @param baseline A length-1 character vector of either "none" (the default),
#' "square", or "trapezoid".
#'
#' @return A length-1 numeric value representing the area under the curve
#' @export
#'
#' @examples
#' trapz(1:10, 1:10)
#' trapz(1:10, 10:1)
#'
#' trapz(1:10, 11:20)
#' trapz(1:10, 11:20, baseline="square")
#' trapz(1:10, 11:20, baseline="trapezoid")
#'
#' x_vals <- seq(-2, 2, length.out=100)
#' trapz(x_vals, dnorm(x_vals))
#'
#' x_vals <- seq(0, pi/2, length.out=100)
#' trapz(x_vals, cos(x_vals))
trapz <- function(rts, ints, baseline="none"){
  # Checks for argument accuracy
  if(!baseline%in%c("none", "square", "trapezoid")){
    stop("Argument 'baseline' must be one of 'none', 'square', or 'trapezoid'")
  }
  if(length(rts)!=length(ints))stop("rts length must be equal to ints length")
  if(!is.numeric(rts))stop("rts vector must be numeric")
  if(!is.numeric(ints))stop("ints vector must be numeric")
  if(length(rts)<2)return(0)

  # Warn if duplicated RT entries (maybe better as any(duplicated(rts)) ?)
  if(length(rts)!=length(unique(rts)))warning("RTs appear to be duplicated")
  # Warn if RTs do not increase monotonically
  # Test stolen from https://stackoverflow.com/q/13093912
  if(!all(rts == cummax(rts)))warning("RTs do not appear to be ordered")
  # Warn with very large peak widths
  if(length(rts)>1000)warning("Most peak widths are less than 1000 scans wide")
  # Warn with negative intensity values
  if(any(ints<0))warning("Detected intensities less than zero")

  # Calculate the area of each trapezoid (area under curve)
  rt_wedges <- rts[-1]-rts[-length(rts)]
  int_wedges <- (ints[-1]+ints[-length(ints)])/2
  auc <- sum(rt_wedges*int_wedges)

  # Calculate area under the baseline (area under baseline)
  if(baseline=="square"){
    aub <- min(ints)*(max(rts)-min(rts))
  } else if(baseline=="trapezoid"){
    minmax_rts <- c(rts[1], rts[length(rts)])
    minmax_ints <- c(ints[1], ints[length(ints)])
    aub <- trapz(minmax_rts, minmax_ints, baseline="none")
  } else {
    aub <- 0
  }

  return(auc-aub)
}



#' Quick plot for MS data
#'
#' Syntactic sugar for a common chromatogram plot. Will use `ggplot2` if
#' available but has a base plot implementation for use even in ultra
#' lightweight situations. Accepts the default MS1 output from `grabMSdata`
#' of a data.table (or base data.frame) with columns for rt (retention time)
#' and int (intensity), usually with filename. Creates a plot of intensity vs
#' retention time with one trace per file (and per level of `color_col`, if
#' that column is discrete). Facets and other customizations can be added to
#' the `ggplot2` output with the usual + syntax.
#'
#' @param MS1_df A data.table with columns named rt and int, and usually filename
#' @param color_col The name of the column to color by. Must be quoted.
#' Discrete columns (character, factor, logical, integer) get one color per
#' level and double columns a continuous gradient (base graphics treat every
#' column as discrete). Defaults to "filename" if that column exists; use NULL
#' for no coloring.
#' @param force_base Boolean option to force base R graphics instead of `ggplot`
#' even if the `ggplot2` package is installed.
#' @param slice_max Boolean option to keep only the most intense point in each
#' scan of each file (and each level of `color_col`), like
#' `dplyr::slice_max(int, by=c(rt, filename, color_col), with_ties=FALSE)`.
#' Useful when a wide m/z window picks up several points per scan, which
#' otherwise draws as a zigzag line.
#'
#' @return If `ggplot2` is installed, a `ggplot` object that can be further
#' modified via additional + commands. Otherwise, NULL (invisibly) and the plot
#' appears via base graphics at the active device.
#' @export
#'
#' @examples
#' \dontshow{data.table::setDTthreads(2)}
#' sample_dir <- system.file("extdata", package = "RaMS")
#' sample_files <- list.files(sample_dir, full.names=TRUE)
#' msdata <- grabMSdata(sample_files[c(3, 5, 6, 7)], grab_what=c("MS1", "BPC"),
#'                      incl_polarity=TRUE)
#'
#' # Base peak chromatogram of every file, all in one color. S30657 switches
#' # polarity every scan, so keep only the positive mode scans to avoid a zigzag
#' qplotMS1data(msdata$BPC[polarity==1], color_col=NULL)
#'
#' # Extracted ion chromatogram of glycine betaine, colored by file
#' qplotMS1data(msdata$MS1[mz%between%pmppm(118.0865, 10)])
#'
#' # A wide m/z window picks up other ions in some scans, drawn as spikes;
#' # slice_max keeps only the biggest point in each scan
#' qplotMS1data(msdata$MS1[mz%between%pmppm(118.0865, 1000)])
#' qplotMS1data(msdata$MS1[mz%between%pmppm(118.0865, 1000)], slice_max=TRUE)
#'
#' # Base graphics work even without ggplot2 installed
#' qplotMS1data(msdata$MS1[mz%between%pmppm(118.0865, 10)], force_base=TRUE)
#'
#' @examplesIf requireNamespace("ggplot2", quietly=TRUE)
#' # The ggplot2 output can be customized as usual, e.g. to mark the
#' # retention time windows used to split the MS2 spectra in ?qplotMS2data
#' library(ggplot2)
#' qplotMS1data(msdata$MS1[mz%between%pmppm(118.0865, 10)]) +
#'   geom_vline(xintercept=c(8.2, 9.5), linetype="dashed")
qplotMS1data <- function(MS1_df, color_col="filename", force_base=FALSE, slice_max=FALSE){
  if(missing(color_col) && !"filename"%in%names(MS1_df))color_col <- NULL
  if(!is.null(color_col)){
    if(!color_col%in%names(MS1_df))stop("Column '", color_col, "' not found in MS1_df")
    # Treat integers (e.g. mz_group IDs) as discrete rather than a gradient
    if(is.integer(MS1_df[[color_col]]))MS1_df[[color_col]] <- factor(MS1_df[[color_col]])
  }
  if(slice_max)MS1_df <- MS1_df[slice_max_int(MS1_df, color_col),]
  if(!force_base && requireNamespace("ggplot2", quietly=TRUE)){
    ggplotMSdata(MS1_df, color_col)
  } else {
    baseplotMSdata(MS1_df, color_col)
  }
}
ggplotMSdata <- function(MS1_df, color_col){
  out_plot <- ggplot2::ggplot(MS1_df) +
    ggplot2::aes(x=rt, y=int) +
    ggplot2::geom_line() +
    ggplot2::labs(x="Retention time (minutes)", y="Intensity")
  trace_syms <- lapply(trace_cols(MS1_df, color_col), as.name)
  if(length(trace_syms)>0){
    out_plot <- out_plot + ggplot2::aes(group=interaction(!!!trace_syms))
  }
  if(!is.null(color_col)){
    out_plot <- out_plot + ggplot2::aes(color=!!as.name(color_col))
  }
  int <- NULL #To prevent R CMD check "notes"  when using ggplot2 syntax
  rt <- NULL
  filename <- NULL
  out_plot
}
baseplotMSdata <- function(MS1_df, color_col){
  if(all(is.na(MS1_df$int)))stop("No intensity data to plot")
  graphics::plot(MS1_df$rt, MS1_df$int, type="n",
                 xlab = "Retention time (minutes)", ylab="Intensity")
  MS1_df <- MS1_df[order(MS1_df$rt),]
  color_map <- base_colors(MS1_df, color_col)
  trace_vals <- lapply(trace_cols(MS1_df, color_col), function(col)MS1_df[[col]])
  traces <- if(length(trace_vals)>0){
    split(seq_len(nrow(MS1_df)), trace_vals, drop=TRUE)
  } else {
    list(seq_len(nrow(MS1_df)))
  }
  for(idx in traces){
    # Each line takes the color of its first point, which is exact for discrete
    # and per-file columns
    graphics::lines(MS1_df$rt[idx], MS1_df$int[idx], col=color_map$colors[idx[1]])
  }
  if(!is.null(color_col)){
    graphics::legend("topright", legend=color_map$labels, col=color_map$legend_cols,
                     lty=1, lwd=2, title=color_col, bty="n", cex=0.8)
  }
  invisible(NULL)
}
# Columns that split the data into separate traces: filename (if present) and
# a discrete color column, so e.g. each mz_group in each file gets its own line
trace_cols <- function(df, color_col){
  unique(c(if("filename"%in%names(df)) "filename",
           if(!is.null(color_col) && !is.numeric(df[[color_col]])) color_col))
}
# Like dplyr::slice_max(int, by=c(rt, <trace_cols>), with_ties=FALSE): flags
# the most intense point in each scan of each trace. order() is stable, so a
# tie goes to whichever point came first
slice_max_int <- function(MS1_df, color_col=NULL){
  group_cols <- c("rt", trace_cols(MS1_df, color_col))
  group_key <- do.call(paste, lapply(group_cols, function(col)MS1_df[[col]]))
  int_ord <- order(MS1_df$int, decreasing=TRUE)
  seq_len(nrow(MS1_df)) %in% int_ord[!duplicated(group_key[int_ord])]
}
# One color per row for the base plots, treating every color column as
# discrete, plus the labels and colors for a legend
base_colors <- function(df, color_col){
  if(is.null(color_col))return(list(colors=rep("black", nrow(df))))
  color_vals <- as.factor(df[[color_col]])
  level_cols <- grDevices::hcl.colors(nlevels(color_vals), palette="Dark 3")
  list(colors=level_cols[color_vals], labels=levels(color_vals), legend_cols=level_cols)
}

#' Quick plot for MS2 data
#'
#' Syntactic sugar for a common fragmentation spectrum plot, drawing each
#' fragment as a vertical line from zero up to its intensity with a point on
#' top. Like \code{\link{qplotMS1data}}, uses `ggplot2` if available and base
#' graphics otherwise. Accepts the default MS2 output from `grabMSdata`, or
#' any data.table or data.frame with columns named fragmz and int. Usually
#' filtered to a single precursor first (e.g. with `premz` and `pmppm`).
#'
#' @param MS2_df A data.table with at least two columns named fragmz and int
#' @param color_col The name of the column to color by. Must be quoted.
#' Discrete columns (character, factor, logical, integer) get one color per
#' level and double columns a continuous gradient (base graphics treat every
#' column as discrete). Defaults to "filename" if that column exists; use NULL
#' for no coloring.
#' @param force_base Boolean option to force base R graphics instead of `ggplot`
#' even if the `ggplot2` package is installed.
#' @param relative Boolean option to scale each spectrum so that its most
#' intense fragment is 100\%, making spectra of very different intensities
#' comparable by shape. A spectrum is all the rows sharing a filename and rt
#' (i.e. one scan), using whichever of those columns are present.
#'
#' @return If `ggplot2` is installed, a `ggplot` object that can be further
#' modified via additional + commands. Otherwise, NULL (invisibly) and the plot
#' appears via base graphics at the active device.
#' @export
#'
#' @examples
#' \dontshow{data.table::setDTthreads(2)}
#' sample_dir <- system.file("extdata", package = "RaMS")
#' sample_files <- list.files(sample_dir, full.names=TRUE)
#' msdata <- grabMSdata(sample_files[7], grab_what="MS2")
#'
#' # All the MS2 scans with glycine betaine's m/z as the precursor
#' betaine_ms2 <- msdata$MS2[premz%between%pmppm(118.0865, 10)]
#' qplotMS2data(betaine_ms2)
#'
#' # Scale each scan to its biggest fragment so smaller scans are visible
#' qplotMS2data(betaine_ms2, relative=TRUE)
#' qplotMS2data(betaine_ms2, relative=TRUE, force_base=TRUE)
#'
#' @examplesIf requireNamespace("ggplot2", quietly=TRUE) && requireNamespace("dplyr", quietly=TRUE)
#' # The three scans come from three different chromatographic peaks (see the
#' # last example in ?qplotMS1data), so split them up and compare their spectra
#' library(ggplot2)
#' library(dplyr)
#' betaine_ms2 %>%
#'   mutate(rt_bin=cut(rt, breaks=c(4, 8.2, 9.5, 15))) %>%
#'   qplotMS2data(relative=TRUE) +
#'   facet_wrap(~rt_bin, ncol=1)
qplotMS2data <- function(MS2_df, color_col="filename", force_base=FALSE, relative=FALSE){
  if(missing(color_col) && !"filename"%in%names(MS2_df))color_col <- NULL
  if(!is.null(color_col)){
    if(!color_col%in%names(MS2_df))stop("Column '", color_col, "' not found in MS2_df")
    # Treat integers (e.g. voltage) as discrete rather than a gradient
    if(is.integer(MS2_df[[color_col]]))MS2_df[[color_col]] <- factor(MS2_df[[color_col]])
  }
  y_lab <- "Intensity"
  if(relative){
    # Percent of the most intense fragment in each spectrum (filename + rt)
    spectrum_vals <- lapply(intersect(c("filename", "rt"), names(MS2_df)), function(col)MS2_df[[col]])
    spectrum_max <- do.call(stats::ave, c(list(MS2_df$int), spectrum_vals, list(FUN=max)))
    MS2_df$int <- MS2_df$int/spectrum_max*100
    y_lab <- "Relative intensity (%)"
  }
  if(!force_base && requireNamespace("ggplot2", quietly=TRUE)){
    ggplotMS2data(MS2_df, color_col, y_lab)
  } else {
    baseplotMS2data(MS2_df, color_col, y_lab)
  }
}
ggplotMS2data <- function(MS2_df, color_col, y_lab){
  out_plot <- ggplot2::ggplot(MS2_df) +
    ggplot2::aes(x=fragmz, y=int, xend=fragmz, yend=0) +
    ggplot2::geom_segment() +
    ggplot2::geom_point() +
    ggplot2::labs(x="Fragment m/z", y=y_lab)
  if(!is.null(color_col)){
    out_plot <- out_plot + ggplot2::aes(color=!!as.name(color_col))
  }
  int <- NULL #To prevent R CMD check "notes"  when using ggplot2 syntax
  fragmz <- NULL
  out_plot
}
baseplotMS2data <- function(MS2_df, color_col, y_lab){
  if(all(is.na(MS2_df$int)))stop("No intensity data to plot")
  graphics::plot(MS2_df$fragmz, MS2_df$int, type="n",
                 ylim=c(0, max(MS2_df$int, na.rm=TRUE)),
                 xlab="Fragment m/z", ylab=y_lab)
  color_map <- base_colors(MS2_df, color_col)
  graphics::segments(x0=MS2_df$fragmz, y0=0, x1=MS2_df$fragmz, y1=MS2_df$int,
                     col=color_map$colors)
  graphics::points(MS2_df$fragmz, MS2_df$int, pch=19, cex=0.7, col=color_map$colors)
  if(!is.null(color_col)){
    graphics::legend("topright", legend=color_map$labels, col=color_map$legend_cols,
                     lty=1, lwd=2, title=color_col, bty="n", cex=0.8)
  }
  invisible(NULL)
}

#' Group m/z values into bins of a specified ppm width
#'
#' This function bins m/z values based on their proximity to each other in m/z
#' space. The algorithm takes the first value in the m/z vector and uses that
#' as the center of a window with a ppm value provided by the user and assigns
#' all m/z values within that window to the same group, then removes those
#' values from consideration and repeats the process until there are no points
#' left to group. This is often used to construct chromatograms from raw MS
#' data that can then be visualized or peakpicked. The function can also drop
#' groups of m/z values if there's not enough points within them or produce
#' only a certain number of groups. Because the algorithm uses the first value
#' in the m/z vector as the window center, it's often a good idea to first
#' sort the values by decreasing intensity.
#'
#' Internally the values are sorted once and each window is located by binary
#' search, so run time scales roughly as n log n rather than n times the number
#' of groups, making it practical on full MS1 tables as well as EIC slices.
#'
#' @param mz_vals A numeric vector of m/z values. NAs are allowed and are
#' returned as NA.
#' @param ppm A length-1 positive numeric specifying the window half-width in
#' ppm, as in \code{pmppm()}.
#' @param min_group_size A length-1 numeric specifying the minimum number of
#' points that must fall within an m/z window for it to be kept as a group.
#' Points in windows smaller than this are returned as NA.
#' @param max_groups A length-1 numeric specifying the maximum number of groups
#' to return. Groups are kept in the order they are found (i.e. by the order
#' of \code{mz_vals}), and windows dropped by \code{min_group_size} do not count
#' toward this limit. If fewer groups exist, all of them are returned.
#' @param mz_factor Logical. If FALSE (the default), groups are identified by
#' consecutive integers in the order they were found. If TRUE, a factor is
#' returned whose levels are the mean m/z of each group's members formatted to
#' 5 decimal places (e.g. "118.08626"), ordered by increasing m/z. Use
#' \code{as.numeric(as.character(x))} to recover the masses, since
#' \code{as.numeric(x)} returns the level codes.
#'
#' @return A vector of the same length as mz_vals specifying the group into
#' which each m/z value was binned: integer if \code{mz_factor = FALSE},
#' factor if \code{mz_factor = TRUE}. Values not assigned to a group are returned as NA.
#' @export
#'
#' @examples
#' \dontshow{data.table::setDTthreads(2)}
#' example_mz_vals <- c(118.0, 118.1, 138.0, 152.0, 118.2, 138.1, 118.1)
#' mz_group(example_mz_vals, ppm = 1)
#' mz_group(example_mz_vals, ppm = 1000)
#' mz_group(example_mz_vals, ppm = 200000)
#'
#' mz_group(example_mz_vals, ppm = 1000, min_group_size = 2)
#' mz_group(example_mz_vals, ppm = 1000, max_groups = 2)
#' mz_group(example_mz_vals, ppm = 1000, mz_factor = TRUE)
#'
#' \dontrun{
#' sample_dir <- system.file("extdata", package = "RaMS")
#' sample_files <- list.files(sample_dir, full.names=TRUE)
#' msdata <- grabMSdata(sample_files[c(3, 5, 6)], grab_what="MS1")
#'
#' grouped_MS1 <- msdata$MS1[mz%between%pmppm(119.0865, 100)][
#'  order(int, decreasing = TRUE)][
#'    ,mz_group:=mz_group(mz, ppm = 5)][]
#' print(grouped_MS1)
#'
#' library(ggplot2)
#' library(dplyr)
#' msdata$MS1[mz%between%pmppm(119.0865, 100)] %>%
#'   arrange(desc(int)) %>%
#'   mutate(mz_group=mz_group(mz, ppm=10)) %>%
#'   ggplot() +
#'   geom_point(aes(x=rt, y=mz, color=factor(mz_group)))
#'
#' msdata$MS1[mz%between%pmppm(119.0865, 100)] %>%
#'   arrange(desc(int)) %>%
#'   mutate(mz_group=mz_group(mz, ppm=5, mz_factor=TRUE)) %>%
#'   qplotMS1data(facet_col = "mz_group")
#' msdata$MS1[mz%between%pmppm(119.0865, 100)] %>%
#'   arrange(desc(int)) %>%
#'   mutate(mz_group=mz_group(mz, ppm=5, max_groups = 2)) %>%
#'   qplotMS1data(facet_col = "mz_group")
#' }
mz_group <- function(mz_vals, ppm, min_group_size=1, max_groups=NULL, mz_factor=FALSE){
  is_num1 <- function(x)is.numeric(x) && length(x)==1 && !is.na(x)
  if(!is.numeric(mz_vals))stop("'mz_vals' must be a numeric vector")
  if(!is_num1(ppm) || ppm<=0)stop("'ppm' must be a length-1 positive numeric")
  if(!is_num1(min_group_size))stop("'min_group_size' must be a length-1 numeric")
  if(!is.null(max_groups)){
    if(!is_num1(max_groups) || max_groups<0)
      stop("'max_groups' must be a length-1 non-negative numeric")
  }
  if(!(isTRUE(mz_factor) || isFALSE(mz_factor)))stop("'mz_factor' must be TRUE or FALSE")

  n <- length(mz_vals)
  mz_groups <- rep(NA_integer_, n)
  has_mz <- which(!is.na(mz_vals))
  max_groups <- if(is.null(max_groups)) Inf else max_groups

  if(length(has_mz) && max_groups>=1){
    # Sort once and find every point's window [mz-ppm, mz+ppm] (inclusive, same
    # arithmetic as pmppm) by binary search, so each pass of the loop below
    # only touches the points inside its own window. The window bounds are
    # looked up in sorted order because findInterval() starts each search from
    # the previous result, which is much faster than random-order lookups
    mz_ord <- has_mz[order(mz_vals[has_mz])]
    sorted_mz <- mz_vals[mz_ord]
    sorted_pos <- integer(n)
    sorted_pos[mz_ord] <- seq_along(mz_ord)
    win_lo <- (findInterval(sorted_mz*(1-ppm/1e6), sorted_mz, left.open=TRUE) + 1L)[sorted_pos]
    win_hi <- findInterval(sorted_mz*(1+ppm/1e6), sorted_mz)[sorted_pos]

    taken <- logical(length(mz_ord))
    group_idx <- 0L
    for(i in has_mz){
      if(taken[sorted_pos[i]])next
      win_idxs <- win_lo[i]:win_hi[i]
      win_idxs <- win_idxs[!taken[win_idxs]]
      taken[win_idxs] <- TRUE
      if(length(win_idxs)>=min_group_size){
        group_idx <- group_idx + 1L
        mz_groups[mz_ord[win_idxs]] <- group_idx
        if(group_idx>=max_groups)break
      }
    }
  }

  if(mz_factor){
    in_group <- !is.na(mz_groups)
    group_means <- rowsum(mz_vals[in_group], mz_groups[in_group], reorder=TRUE)[,1]/
      tabulate(mz_groups[in_group])
    # Build the factor from integer codes directly; factor() would coerce
    # every value to character first, which is slow on long vectors. Groups
    # whose means round to the same label share a level
    mz_labels <- sprintf("%.5f", group_means)
    mz_levels <- unique(mz_labels[order(group_means)])
    mz_groups <- structure(match(mz_labels, mz_levels)[mz_groups],
                           levels=mz_levels, class="factor")
  }

  mz_groups
}



#' Edit mzML/mzXML file retention times
#'
#' This function uses the basic XML parsing of RaMS to modify the retention
#' times of MS scans within the mzML/mzXML files. This method can be useful for
#' performing RT correction using one platform and then peakpicking and
#' correspondence on another. The basic method is simply replacing the scan's
#' recorded retention time value with an arbitrary one of your choosing. This
#' function is vectorized to handle multiple files, while the internal
#' `editMzmlRTs()` and `editMzxmlRTs()` do most of the heavy lifting. Note that
#' the seconds vs minutes must be closely monitored here - the unit should be
#' the same as the one in the file originally.
#'
#' @param files Vector of filenames (including the relative/absolute path)
#' @param new_rt_list Nested of new retention times. One entry in the list for
#' each file (in the same order as the files), each containing a vector of new
#' retention times. RT vectors can be equal to either every scan or just every
#' MS1 scan. If only the MS1 scans are provided in a file with additional MS
#' levels, MSn scan RTs will be interpolated according to interp_method (below)
#' @param new_filenames Vector of filenames (including relative/absolute paths)
#' describing where the edited files should be written out. Can be the same
#' as files but will throw a warning and append _rtcor to each file unless
#' `overwrite = TRUE` (below)
#' @param interp_method Either "linear" or "constant". Describes the way that
#' MSn retention times should be handled when only the MS1 values are provided.
#' "linear" (the default) means that the spacing will be preserved, while
#' "constant" will use the associated MS1 scan RT for all MSn scans, allowing
#' an easy method of linking the MSn to the MS1.
#' @param overwrite Boolean. Controls whether files are overwritten in place
#' if `new_filenames` is not provided.
#'
#' @return Invisibly, the names of the edited files.
#' @export
#'
#' @examples
#' \dontrun{
#' # Setup (allows running on CRAN computers)
#' example_dir <- tempdir()
#' rams_dir <- system.file("extdata", package = "RaMS")
#' file.copy(list.files(rams_dir, pattern = "LB.*mzML", full.names = TRUE), example_dir)
#' mzMLs <- list.files(example_dir, pattern = "LB.*mzML", full.names = TRUE)
#'
#' library(xcms)
#' library(RaMS)
#'
#' register(SerialParam())
#' xcms_obj <- readMSData(mzMLs, msLevel. = 1, mode = "onDisk")
#' cwp <- CentWaveParam(ppm = 5, peakwidth = c(20, 80))
#' xcms_peakpicked <- findChromPeaks(xcms_obj, param = cwp)
#' xcms_rtcor <- adjustRtime(xcms_peakpicked, param = ObiwarpParam())
#'
#' # Extract the adjusted RTs from the XCMS object
#' new_rts <- split(rtime(xcms_rtcor)/60, fromFile(xcms_rtcor))
#' # Apply the retention time correction to the new files
#' mzMLs_rtcor <- editMSfileRTs(mzMLs, new_rt_list = new_rts)
#'
#' # Contrast the two chromatograms to see the peaks aligned
#' qplotMS1data(grabMSdata(mzMLs)$MS1[mz%between%pmppm(104.1073, 10)])
#' qplotMS1data(grabMSdata(mzMLs_rtcor)$MS1[mz%between%pmppm(104.1073, 10)])
#'
#' # Cleanup
#' file.remove(mzMLs)
#' file.remove(mzMLs_rtcor)
#' }
editMSfileRTs <- function(files, new_rt_list, new_filenames=NULL,
                          interp_method="linear", overwrite=FALSE){
  existing_files <- file.exists(files)
  if(!all(existing_files)){
    stop_msg <- paste0(
      "Not all files found, missing\n",
      paste(files[!existing_files], collapse = "\n")
    )
    stop(stop_msg)
  }
  if(length(files)!=length(new_rt_list)){
    stop("'files' and 'new_rt_list' must be the same length")
  }
  if(is.null(new_filenames))new_filenames <- files
  if(any(new_filenames%in%files) & overwrite==FALSE){
    warn_msg <- paste0(
      "New files would overwrite existing files with overwrite=FALSE\n",
      "Appending _rtcor to each"
    )
    warning(warn_msg)
    new_filenames <- gsub("\\.(?=mzX?ML)", replacement = "_rtcor\\.",
                          new_filenames, perl = TRUE)
  }
  if(length(files)!=length(new_rt_list)){
    stop("'files' and 'new_rt_list' must be the same length")
  }

  mzML_file_idxs <- grep("\\.mzML(\\.gz)?$", files)
  mapply(editMzmlRTs, files[mzML_file_idxs], new_filenames[mzML_file_idxs],
         new_rt_list[mzML_file_idxs], interp_method)

  mzXML_file_idxs <- grep("\\.mzXML(\\.gz)?$", files)
  mapply(editMzxmlRTs, files[mzXML_file_idxs], new_filenames[mzXML_file_idxs],
         new_rt_list[mzXML_file_idxs], interp_method)

  invisible(new_filenames)
}
editMzmlRTs <- function(filename, output_filename, new_rts, interp_method="linear"){
  xml_data <- xml2::read_xml(filename)
  checkFileType(xml_data, "mzML")
  scan_xpath <- '//d1:spectrum[d1:cvParam[@name="ms level"]]'
  scan_nodes <- xml2::xml_find_all(xml_data, scan_xpath)

  if(length(scan_nodes)==length(new_rts)){
    new_rt_vals <- new_rts
  } else {
    ms1_xpath <- '//d1:spectrum[d1:cvParam[@name="ms level" and @value="1"]]'
    ms1_nodes <- xml2::xml_find_all(xml_data, ms1_xpath)

    if(length(ms1_nodes)!=length(new_rts)){
      stop_msg <- paste0(
        "new_rts must be equal in length to either the total number of scans (",
        length(scan_nodes), ") or the number of MS1 scans (", length(ms1_nodes),
        ") but does not seem to match either (length = ", length(new_rts), ")"
      )
      stop(stop_msg)
    }

    rt_xpath <- 'd1:scanList/d1:scan/d1:cvParam[@name="scan start time"]'
    all_rt_nodes <- xml2::xml_find_all(scan_nodes, rt_xpath)
    all_init_rts <- as.numeric(xml2::xml_attr(all_rt_nodes, "value"))
    ms1_rt_nodes <- xml2::xml_find_all(ms1_nodes, rt_xpath)
    ms1_init_rts <- as.numeric(xml2::xml_attr(ms1_rt_nodes, "value"))

    new_rt_vals <- approx(x = ms1_init_rts, y = new_rts, xout = all_init_rts,
                          method=interp_method)$y
  }

  rt_xpath <- 'd1:scanList/d1:scan/d1:cvParam[@name="scan start time"]'
  rt_nodes <- xml2::xml_find_all(scan_nodes, rt_xpath)
  xml2::xml_attr(rt_nodes, "value") <- new_rt_vals

  # Add note that RaMS was used to edit RTs (stolen from minifyMSFunctions.R)
  proclist_node <- xml2::xml_find_all(xml_data, "//d1:dataProcessingList")
  xml2::xml_add_child(proclist_node, "dataProcessing", id="RaMS_R_package")
  proc_node <- xml2::xml_find_all(proclist_node, '//dataProcessing[@id="RaMS_R_package"]')
  xml2::xml_add_child(proc_node, "processingMethod", order=0, softwareRef="RaMS")
  meth_node <- xml2::xml_find_all(proclist_node, '//processingMethod[@order="0"]')
  xml2::xml_add_child(meth_node, "userParam", cvRef="MS", accession="MS:1009003",
                      name="RT editing via RaMS", value="")
  process_count <- as.numeric(xml2::xml_attr(proclist_node, "count"))+1
  xml2::xml_attr(proclist_node, "count") <- process_count
  softlist_node <- xml2::xml_find_all(xml_data, "//d1:softwareList")
  xml2::xml_add_child(softlist_node, "software", id="RaMS",
                      version=as.character(packageVersion("RaMS")))
  soft_node <- xml2::xml_find_all(softlist_node, 'software[@id="RaMS"]')
  xml2::xml_add_child(soft_node, "userParam", name="RaMS R package")
  software_count <- as.numeric(xml2::xml_attr(softlist_node, "count"))+1
  xml2::xml_attr(softlist_node, "count") <- software_count

  # Write the file out
  xml2::write_xml(xml_data, file = output_filename)
}
editMzxmlRTs <- function(filename, output_filename, new_rts, interp_method="linear"){
  xml_data <- xml2::read_xml(filename)
  checkFileType(xml_data, "mzXML")
  scan_xpath <- '//d1:scan'
  scan_nodes <- xml2::xml_find_all(xml_data, scan_xpath)

  if(length(scan_nodes)==length(new_rts)){
    new_rt_vals <- new_rts
  } else {
    ms1_xpath <- '//d1:scan[@msLevel="1" and @peaksCount>0]'
    ms1_nodes <- xml2::xml_find_all(xml_data, ms1_xpath)

    if(length(ms1_nodes)!=length(new_rts)){
      stop_msg <- paste0(
        "new_rts must be equal in length to either the total number of scans (",
        length(scan_nodes), ") or the number of MS1 scans (", length(ms1_nodes),
        ") but does not seem to match either (length = ", length(new_rts), ")"
      )
      stop(stop_msg)
    }

    all_rt_attrs <- xml2::xml_attr(scan_nodes, "retentionTime")
    all_init_prefixes <- gsub("[0-9].*", "", all_rt_attrs)
    all_init_suffixes <- gsub(".*[0-9]", "", all_rt_attrs)
    all_init_rts <- as.numeric(gsub("PT|S", "", all_rt_attrs))
    ms1_rt_attrs <- xml2::xml_attr(ms1_nodes, "retentionTime")
    ms1_init_rts <- as.numeric(gsub("PT|S", "", ms1_rt_attrs))

    new_rt_vals <- approx(x = ms1_init_rts, y = new_rts, xout = all_init_rts,
                          method=interp_method)$y
    new_rt_vals <- paste0(all_init_prefixes, new_rt_vals, all_init_suffixes)
  }

  rt_nodes <- xml2::xml_find_all(scan_nodes, scan_xpath)
  xml2::xml_attr(rt_nodes, "retentionTime") <- new_rt_vals

  # Add note that RaMS was used to edit RTs (stolen from minifyMSFunctions.R)
  proclist_node <- xml2::xml_find_all(xml_data, "//d1:dataProcessing")
  xml2::xml_add_sibling(proclist_node, "dataProcessing")
  proc_node <- xml2::xml_find_first(xml_data, "//dataProcessing")
  xml2::xml_add_child(proc_node, "software", type="RT editing", name="RaMS",
                      version=as.character(packageVersion("RaMS")))

  # Write the file out
  xml2::write_xml(xml_data, file = output_filename)
}
