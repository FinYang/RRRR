% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RRR.R
\name{RRR}
\alias{RRR}
\title{Reduced-Rank Regression using Gaussian MLE}
\usage{
RRR(y, x, z = NULL, mu = TRUE, r = 1)
}
\arguments{
\item{y}{Matrix of dimension N*P. The matrix for the response variables. See \code{Detail}.}

\item{x}{Matrix of dimension N*Q. The matrix for the explanatory variables to be projected. See \code{Detail}.}

\item{z}{Matrix of dimension N*R. The matrix for the explanatory variables not to be projected. See \code{Detail}.}

\item{mu}{Logical. Indicating if a constant term is included.}

\item{r}{Integer. The rank for the reduced-rank matrix \eqn{AB'}. See \code{Detail}.}
}
\value{
A list of the estimated parameters of class \code{RRR}.
\describe{
\item{spec}{The input specifications. \eqn{N} is the sample size.}
\item{mu}{The estimated constant vector. Can be \code{NULL}.}
\item{A}{The estimated exposure matrix.}
\item{B}{The estimated factor matrix.}
\item{D}{The estimated coefficient matrix of \code{z}. Can be \code{NULL}.}
\item{Sigma}{The estimated covariance matrix of the innovation distribution.}
}
}
\description{
Gaussian Maximum Likelihood Estimation method for Reduced-Rank Regression.
This method is not robust in the sense that it assumes a Gaussian distribution
for the innovations
which does not take into account the heavy-tailedness of the true distribution and
outliers.
}
\details{
The formulation of the reduced-rank regression is as follow:
\deqn{y = \mu +AB'  x + D z+innov,}
where for each realization \eqn{y} is a vector of dimension \eqn{P} for the \eqn{P} response variables,
\eqn{x} is a vector of dimension \eqn{Q} for the \eqn{Q} explanatory variables that will be projected to
reduce the rank,
\eqn{z} is a vector of dimension \eqn{R} for the \eqn{R} explanatory variables
that will not be projected,
\eqn{\mu} is the constant vector of dimension \eqn{P},
\eqn{innov} is the innovation vector of dimension \eqn{P},
\eqn{D} is a coefficient matrix for \eqn{z} with dimension \eqn{P*R},
\eqn{A} is the so called exposure matrix with dimension \eqn{P*r}, and
\eqn{B} is the so called factor matrix with dimension \eqn{Q*r}.
The matrix resulted from \eqn{AB'} will be a reduced rank coefficient matrix with rank of \eqn{r}.
The function estimates parameters \eqn{\mu}, \eqn{A}, \eqn{B}, \eqn{D}, and \eqn{Sigma}, the covariance matrix of
the innovation's distribution, assuming the innovation has a Gaussian distribution.
}
\examples{
set.seed(2222)
data <- RRR_sim()
res <- RRR(y=data$y, x=data$x, z = data$z)
res
}
\references{
S. Johansen, "Estimation and Hypothesis Testing of Cointegration Vectors in Gaussian Vector Autoregressive Models,"Econometrica, vol. 59,p. 1551, Nov. 1991.
}
\seealso{
For robust reduced-rank regression estimation see function \code{\link{RRRR}}.
}
\author{
Yangzhuoran Yang
}
