
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
#' @export
#'
#' @examples
#' pmppm(100, 5)
#' pmppm(1000000, 5)
#' pmppm(118.0865, 2.5)
#' pmppm(892.535313, 10)
pmppm <- function(mass, ppm=4){
  if(length(mass)>1)stop("Vectorized masses are not supported in pmppm")
  if(length(ppm)>1)stop("Vectorized ppms are not supported in pmppm")
  if(!is.numeric(mass))stop("'mass' must be numeric")
  if(!is.numeric(ppm))stop("'ppm' must be numeric")
  if(ppm<0)warning("'ppm' should be positive")
  c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))
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
#' and int (intensity) as well as filename. Creates a plot of intensity vs
#' retention time with one trace per file. A few additional `ggplot2` arguments
#' are also made available for easy coloring or facetting by providing the
#' name of the associated column to the `color_col` and `facet_col` arguments,
#' respectively.
#'
#' @param MS1_df A data.table with at least three columns named rt, int, and filename
#' @param color_col The name of the column to color by. Must be quoted.
#' @param facet_col The name of the column to facet by. Must be quoted.
#' @param facet_args Since the call to facet_wrap is within the function, you
#' can provide additional facet customization arguments here as a list. Although
#' if you're starting to fiddle with facets you'll probably be better served by
#' the proper `ggplot` call.
#' @param force_base Boolean option to force base R graphics instead of `ggplot`
#' even if the `ggplot2` package is installed.
#'
#' @return If `ggplot2` is installed, a `ggplot` object that can be further
#' modified via additional + commands. Otherwise, NULL and the plot appears
#' via base graphics at the active device.
#' @export
#'
#' @examples
#' test_df <- expand.grid(rt=rep(1:100, length.out=1000))
#' test_df$int <- rep(dnorm(seq(-10, 10, length.out=100)), 10)*10+runif(1000)
#' test_df$filename <- rep(LETTERS[1:10], each=100)
#' qplotMS1data(test_df)
#'
#' test_df$startime <- rep(gl(2, 5, labels = c("Morn", "Eve")), each=100)
#' qplotMS1data(test_df, color_col="startime", facet_col="startime")
#' qplotMS1data(test_df, color_col="startime", facet_col="startime",
#'             facet_args=list(ncol=2, scales="free"))
#'
#' # Using data from the `grabMSdata` function:
#' sample_dir <- system.file("extdata", package = "RaMS")
#' sample_files <- list.files(sample_dir, full.names=TRUE)
#' msdata <- grabMSdata(sample_files[c(3, 5, 6)], grab_what="MS1")
#' qplotMS1data(msdata$MS1[mz%between%pmppm(118.0865)])
qplotMS1data <- function(MS1_df, color_col=NULL, facet_col=NULL,
                         facet_args=list(ncol=1), force_base=FALSE){
  if(requireNamespace("ggplot2", quietly=TRUE) & !force_base){
    ggplotMSdata(MS1_df, color_col, facet_col, facet_args)
  } else {
    if(!is.null(color_col))warning("Argument 'color_col' is currently only available via ggplot2")
    if(!is.null(facet_col))warning("Argument 'facet_col' is currently only available via ggplot2")
    baseplotMSdata(MS1_df)
  }
}
ggplotMSdata <- function(MS1_df, color_col, facet_col, facet_args){
  out_plot <- ggplot2::ggplot(MS1_df) +
    ggplot2::aes(x=rt, y=int, group=filename) +
    ggplot2::geom_line() +
    ggplot2::labs(x="Retention time (minutes)", y="Intensity")
  if(!is.null(color_col)){
    out_plot <- out_plot + ggplot2::aes(color=get(color_col)) + ggplot2::labs(color=color_col)
  }
  if(!is.null(facet_col)){
    arg_list <- c(list(facets=ggplot2::vars(get(facet_col))), facet_args)
    out_plot <- out_plot + do.call(what=ggplot2::facet_wrap, args = arg_list)
  }
  out_plot
}
baseplotMSdata <- function(MS1_df){
  plot(MS1_df$rt, MS1_df$int, type="n",
       xlab = "Retention time (minutes)", ylab="Intensity")
  ordered_df <- MS1_df[order(MS1_df$rt),]
  spldf <- split(ordered_df, ordered_df$filename)
  sapply(spldf, function(spldf_i)lines(x=spldf_i$rt, y=spldf_i$int))
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
#' @param mz_vals A numeric vector of m/z values
#' @param ppm A length-1 numeric vector specifying the desired window size in ppm
#' @param min_group_size A length-1 numeric vector specifying the minimum number
#' of points that must fall within an m/z window to be assigned a group number
#' @param max_groups A length-1 numeric vector specifying the maximum number
#' of total groups to assign.
#'
#' @return A numeric vector of the same length as mz_vals specifying the group
#' into which each m/z value was binned. Values not assigned to a group are
#' returned as NAs.
#' @export
#'
#' @examples
#' example_mz_vals <- c(118.0, 118.1, 138.0, 152.0, 118.2, 138.1, 118.1)
#' mz_group(example_mz_vals, ppm = 1)
#' mz_group(example_mz_vals, ppm = 1000)
#' mz_group(example_mz_vals, ppm = 200000)
#'
#' mz_group(example_mz_vals, ppm = 1000, min_group_size = 2)
#' mz_group(example_mz_vals, ppm = 1000, max_groups = 2)
#'
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
#'   mutate(mz_group=mz_group(mz, ppm=5)) %>%
#'   qplotMS1data(facet_col = "mz_group")
#' msdata$MS1[mz%between%pmppm(119.0865, 100)] %>%
#'   arrange(desc(int)) %>%
#'   mutate(mz_group=mz_group(mz, ppm=5, max_groups = 2)) %>%
#'   qplotMS1data(facet_col = "mz_group")
mz_group <- function(mz_vals, ppm, min_group_size=0, max_groups=NULL){
  if(!is.numeric(mz_vals))stop("'mz_vals' must be a numeric vector")
  if(!is.numeric(ppm))stop("'ppm' must be a length-1 numeric")
  if(!is.numeric(min_group_size))stop("'ppm' must be a length-1 numeric")
  if(!is.null(max_groups)){
    if(!is.numeric(max_groups))stop("'max_groups' must be a length-1 numeric")
  }

  mz_groups <- numeric(length(mz_vals))
  needs_group <- !logical(length(mz_vals))

  if(is.null(max_groups)){
    group_idx <- 1L
    while(sum(needs_group)>0){
      mz_idxs <- mz_vals%between%pmppm(mz_vals[needs_group][1], ppm)
      if(sum(mz_idxs[needs_group])>min_group_size){
        mz_groups[mz_idxs&needs_group] <- group_idx
        group_idx <- group_idx + 1L
      }
      needs_group[mz_idxs] <- FALSE
    }
  } else {
    for(group_idx in seq_len(max_groups)){
      mz_idxs <- mz_vals%between%pmppm(mz_vals[needs_group][1], ppm)
      if(sum(mz_idxs[needs_group])>min_group_size){
        mz_groups[mz_idxs&needs_group] <- group_idx
      }
      needs_group[mz_idxs] <- FALSE
    }
  }
  mz_groups[mz_groups==0] <- NA
  mz_groups
}
