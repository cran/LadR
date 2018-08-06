
#main author: Kévin Allan Sales Rodrigues

#### algorithm to make inference in LAD models

#' Fitting LAD Models
#'
#' @param y A vector with response variables.
#' @param x A matrix or vector with explanatory variables.
#' @param intercept TRUE for a model with intercept and FALSE for a model without intercept.
#' @return
#'    list defining the regression (compare with function \code{\link{lsfit}}).
#' \item{coefficients}{vector of coefficients.}
#' \item{residuals}{residuals from the fit.}
#' \item{message}{vector of one or two character strings stating whether a
#'  non-unique solution is possible, or if the x matrix was found to be rank deficient.}
#'
#' @references
#' Barrodale, I., and Roberts, F.D.K. (1973).
#'An improved algorithm for discrete L1 linear approximations.
#'\emph{SIAM Journal of Numerical Analysis} \bold{10}, 839-848.
#'
#'Barrodale, I., and Roberts, F.D.K. (1974).
#'Solution of an overdetermined system of equations in the L1 norm.
#'\emph{Communications of the ACM} \bold{17}, 319-320.
#'
#'Bloomfield, P., and Steiger, W.L. (1983).
#'\emph{Least Absolute Deviations: Theory, Applications, and Algorithms.}
#'Birkhauser, Boston, Mass.
#'
#'
#' @details
#' The Barrodale-Roberts algorithm, which is a specialized linear programming
#' algorithm, is used.
#'
#'
#' @examples
#' ### Using stackloss data
#'
#' ladfit(stack.x, stack.loss, intercept =TRUE)



ladfit = function (x, y, intercept = TRUE) {
tolerance = 1e-07
print.it = TRUE
    xn = if (is.matrix(x))
        dimnames(x)[[2]]
    else "X"
    x = as.matrix(x)
    dx = dim(x)
    vars = dx[2]
    if (length(xn) == 0)
        xn = paste("X", 1:vars, sep = "")
    if (intercept) {
        vars = vars + 1
        xn = c("Intercept", xn)
        x = cbind(1, x)
    }
    obs = dx[1]
    if (obs != length(y))
        stop("Y and X have different number of observations!")
    if (obs <= vars)
        stop("Less observations than variables")
    storage.mode(y) <- "double"
    z <- matrix(0, nrow = obs + 2, ncol = vars + 2)
    z[1:obs, 1:vars] <- x
    storage.mode(z) <- "double"
    fit <- .Fortran("l1", n = as.integer(obs), p = as.integer(vars),
        n2 = as.integer(obs + 2), p2 = as.integer(vars + 2),
        z = z, y = y, tol = as.double(tolerance), cf = double(vars),
        resid = double(obs), work = integer(obs))
    minimum <- fit$z[obs + 1, vars + 1]
    rank = fit$z[obs + 1, vars + 2]
    info = fit$z[obs + 2, vars + 1]
    numIter = fit$z[obs + 2, vars + 2]
    fitted = drop(x %*% fit$cf)
    msg = character(0)
    if (info == 0) {
        msg = "Non-unique solution possible!"
        if (print.it)
            warning(msg)
    }
    else if (info != 1)
        stop("Premature termination!!")
    if (rank != vars) {
        msg = c(msg, paste("Matrix not of full rank, apparent rank",
            rank))
        if (print.it)
            warning(msg[length(msg)])
    }
    out = list(coefficients = fit$cf, minimum = minimum, fitted.values = fitted,
        residuals = fit$resid, rank = rank, numIter = numIter,
        info = info)
    names(out$coefficients) <- xn
    if (length(msg))
        out$message = list(msg)
    out
}
