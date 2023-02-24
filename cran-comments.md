This patch aims to fix CRAN Package Check warning and note.
Also updated url.

## Test environments
* local macOS 13.2.1 (Apple M1 Pro) install, R 4.2.2
* maxOS-latest (release) x86_64-apple-darwin17.0 (64-bit) (GitHub actions)
* windows-latest (release) x86_64-w64-mingw32 (64-bit) (GitHub actions)
* ubuntu-latest (devel) x86_64-pc-linux-gnu (64-bit) (GitHub actions)
* ubuntu-latest (release) x86_64-pc-linux-gnu (64-bit) (GitHub actions)
* ubuntu-latest (oldrel-1) x86_64-pc-linux-gnu (64-bit) (GitHub actions)
*	Windows Server 2022, R-devel, 64 bit (on Rhub)
* R version 4.2.2 (2022-10-31 ucrt) x86_64-w64-mingw32 (64-bit) (win-builder)
* R Under development (unstable) (2023-02-23 r83894 ucrt) x86_64-w64-mingw32 (64-bit) (win-builder)

## R CMD check results

0 errors | 0 warnings | 4 note

* hecking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'
  
Only appears on Windows R-devel version. 
From [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX and can be ignored.

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
  Skipping checking math rendering: package 'V8' unavailable

Only on Fedora Linux, R-devel, clang, gfortran (on Rhub).
Cannot manage dependencies used in check. The check is ok on other platforms. Can be ignored.


* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Yangzhuoran Fin Yang <yangyangzhuoran@gmail.com>'
  New maintainer:
  Yangzhuoran Fin Yang <yangyangzhuoran@gmail.com>
  Old maintainer(s):
  Yangzhuoran Yang <Fin.Yang@monash.edu>

Same maintainer. Changed email to a more sustainable one.

* Found the following (possibly) invalid URLs:
  URL: https://www.jstor.org/stable/2938278?origin=crossref
    From: inst/doc/RRRR.html
    Status: 403
    Message: Forbidden

The URLs are valid.

## Downstream dependencies

There are currently no downstream dependencies for this package.

