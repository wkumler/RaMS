
## Resubmission

This is a resubmission. In this version I have reduced the example runtimes to
below 5s (I hope) for tests and examples and set data.table to use only two
threads to avoid CPU time > 2.5 times elapsed time

## Test environments
* local R installation, R 4.3.2
* windows-latest (release) (via GitHub Actions)
* macOS-latest (release) (via GitHub Actions)
* ubuntu-latest (release) (via GitHub Actions)

## R CMD check results

0 errors | 0 warnings | 0 notes

## revdepcheck results

I checked 1 reverse dependency (from CRAN), and had no problems compiling it 
according to tools::check_packages_in_dir
