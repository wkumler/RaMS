% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grabMzmlFunctions.R
\name{grabMzmlEncodingData}
\alias{grabMzmlEncodingData}
\title{Helper function to extract mzML file encoding data}
\usage{
grabMzmlEncodingData(xml_data)
}
\arguments{
\item{xml_data}{mzML data as parsed by xml2}
}
\value{
A list of values used by other parsing functions, currently
compression, mz_precision, int_precision
}
\description{
Helper function to extract mzML file encoding data
}
