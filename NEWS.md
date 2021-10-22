# RaMS (development version)

# RaMS 1.1.0
## The minification release
 - Added `minifyMSdata`, a function that shrinks MS files by excluding masses
 - Restructured demo files using the new `minify*` functions to include more useful data
 - Added `prefilter` argument to `grabMSdata` that drops masses below a user-provided intensity threshold
 - Added quick-start poster to README.md
 - Added `grabAccessionData`, a function that extracts arbitrary MS metadata by HUPO-PSI accession number
 - Squashed some bugs

# RaMS 1.0.0
## The first CRAN release
 - Enabled rtrange parameter for mzXML files
 - All retention times are now converted to minutes rather than seconds
 - Added several vignettes and switched back to knitr for PDF support, abandoning plotly things
   - RaMS and friends: some examples of moving data between R and other languages
   - Basic integrations: some rudimentary code that can be used to manually integrate peaks
 - Streamlined README.md

# RaMS 0.3.0
## The documentation update
 - Expanded documentation for individual functions, both exposed and unexposed
 - Expanded README.md to include simpler demos and link to vignette
 - Created vignette to demonstrate the nuances of `RaMS`
 - Created tests to stabilize future package development
 - Streamlined example mzML and mzXML files

# RaMS 0.2.0
## The user-friendly update
 - Combined mzML and mzXML handling into one wrapper, `grabMSdata`, which also handles multiple files simultaneously
 - Created README.md and included several demos of package functionality
 - Added several "minified" mzML and mzXML files for demos and examples

# RaMS 0.1.0
## The first update
 - Defined major functions, constructed package skeleton
