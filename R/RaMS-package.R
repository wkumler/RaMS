#' @md
#' @details
#' The main function of RaMS is [grabMSdata()], which handles most use-cases and
#' automatically detects file types and is vectorized over multiple files. If
#' more control over the file reading is desired, both [grabMzmlData()] and
#' [grabMzxmlData()] have been exposed to the user, and [grabAccessionData()] is
#' available if a specific accession value from an mzML file is desired.
#'
#' Other useful functions in the package include [minifyMSdata()] and
#' [tmzmlMaker()] which have their own documentation pages and vignettes. A few
#' small helper functions are also included: [pmppm()] for building a ppm range
#' from a single mass (and its vectorized form [pmppm_vec()]); [trapz()] for
#' peak integration given RT/int vectors, and [qplotMS1data()]
#' and [qplotMS2data()] because I got tired of typing out the same plotting
#' commands. [`%between%`][reexports] and [`%inrange%`][reexports] are imported
#' from `data.table` for ease of use.
#'
#' See the package intro on GitHub at https://github.com/wkumler/RaMS and
#' explore the vignette with \code{vignette("Intro-to-RaMS", package = "RaMS")}
#'
#' @keywords internal
"_PACKAGE"

#' @import xml2
#' @import data.table
#' @import utils
#' @importFrom base64enc base64decode base64encode
#' @importFrom stats approx
NULL
