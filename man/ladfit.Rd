% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ladfit.R
\name{ladfit}
\alias{ladfit}
\title{Fitting LAD Models}
\usage{
ladfit(x, y, intercept = TRUE)
}
\arguments{
\item{x}{A matrix or vector with explanatory variables.}

\item{y}{A vector with response variables.}

\item{intercept}{TRUE for a model with intercept and FALSE for a model without intercept.}
}
\value{
list defining the regression (compare with function \code{\link{lsfit}}).
\item{coefficients}{vector of coefficients.}
\item{residuals}{residuals from the fit.}
\item{message}{vector of one or two character strings stating whether a
 non-unique solution is possible, or if the x matrix was found to be rank deficient.}
}
\description{
Fitting LAD Models
}
\details{
The Barrodale-Roberts algorithm, which is a specialized linear programming
algorithm, is used.
}
\examples{
### Using stackloss data

ladfit(stack.x, stack.loss, intercept =TRUE)
}
\references{
Barrodale, I., and Roberts, F.D.K. (1973).
An improved algorithm for discrete L1 linear approximations.
\emph{SIAM Journal of Numerical Analysis} \bold{10}, 839-848.

Barrodale, I., and Roberts, F.D.K. (1974).
Solution of an overdetermined system of equations in the L1 norm.
\emph{Communications of the ACM} \bold{17}, 319-320.

Bloomfield, P., and Steiger, W.L. (1983).
\emph{Least Absolute Deviations: Theory, Applications, and Algorithms.}
Birkhauser, Boston, Mass.
}
