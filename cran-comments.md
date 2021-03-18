## Test environments
* local R installation, R 4.0.4
* windows-latest (release) (via GitHub Actions)
* macOS-latest (release) (via GitHub Actions)
* ubuntu-latest (release) (via GitHub Actions)

## R CMD check results

0 errors | 0 warnings | 1 note

* **This is a new release.**

On some builds, an additional note is returned:

* **Examples with CPU (user + system) or elapsed time > 5s.** 
The examples for grabMSdata typically take a while because there are a large
number of examples demonstrating multiple parameters and because the function
typically reads a large amount of data from disk. Both of these are expected
behaviors.

## References describing the methods in package
* There are no references describing the methods this package. These are in
progress but existence on CRAN is desired prior to publication of these methods.



