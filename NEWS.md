# RRRR 1.1.0.9000

* Dependency of `expm` changed to `pracma`. The package is used for calculating the square root of a matrix. `expm` is scheduled for archival on 2020-07-29 for unfixed check problem.

# RRRR 1.1.0

* Improved documentation
* Bug fix: update.RRRR now includes previous history
* Improved safety: added more checks in RRRR and ORRRR
* Bug fix: plot.RRRR run time x aes wrong labels are now fixed.
* Bug fix: coef.RRR method now handles different Q and R
* Bug fix: RRR_data class now print different Q and R
* Better print method: extra dimensions are now printed as <unspecified>.
* The package support of the progress bar is now changed from package `dplyr` to package `lazybar`. The progress bar functionality will be deprecated in `dplyr` 1.0.0. `lazybar` provide a better estimated remaining time using drift forecast method.


# RRRR 1.0.0

* First release
