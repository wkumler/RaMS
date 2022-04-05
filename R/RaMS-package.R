#' @md
#' @details
#' The main function of RaMS is [grabMSdata()], which handles most use-cases and
#' automatically detects file types and is vectorized over multiple files. If
#' more control over the file reading is desired, both [grabMzmlData()] and
#' [grabMzxmlData()] have been exposed to the user, and [grabAccessionData()] is
#' available if a specific accession value from an mzML file is desired.
#'
#' Other useful functions in the package include [minifyMSdata()] and
#' [tmzmlMaker()] which have their own documentation pages and vignettes. Two
#' small helper functions are also included, [pmppm()] and [between()] (with the
#' alias `%between%`) imported from `data.table`.
#'
#' See the package intro on GitHub at https://github.com/wkumler/RaMS and
#' explore the vignettes with \code{vignette("help", package = "mypkg")}
#'
#' @keywords internal
"_PACKAGE"
