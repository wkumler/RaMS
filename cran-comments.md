
## Resubmission

This is a resubmission, updating the package from v1.4.0 to v1.4.3 with some
additional bug fixes. I've also tested with the system environment variable
_R_CHECK_DEPENDS_ONLY_ set to true to avoid the problems in my last submission.

## Test environments
* local R installation, R 4.4.0
* windows-latest (release) (via GitHub Actions)
* macOS-latest (release) (via GitHub Actions)
* ubuntu-latest (release) (via GitHub Actions)

## R CMD check results

0 errors | 0 warnings | 0 notes

## revdepcheck results

I checked 1 reverse dependency (from CRAN), and had no problems compiling it 
according to tools::check_packages_in_dir
