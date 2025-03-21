# Flexible Covariate Modelling in Space and Time with the Generalised
# Extreme-value Distribution and Allowing Parameter Constraints
{#
  # Description:
  #
  #      Maximum-likelihood fitting of a Generalised Extreme-value
  #      distribution to block maxima with covariates included in the
  #      location, scale and shape parameters, plus other user-defined
  #      parameter constraints.
  #
  # Usage:
  #
  #      fgev(data, start, fpar, xpar, std.err = TRUE, ...)
  #
  # Arguments:
  #
  #     data: A numeric m-by-n matrix of block maxima.
  #
  #    start: A numeric vector of starting values for the parameters.
  #
  #     fpar: A function taking a parameter vector as its first argument
  #           and `xpar' as its second argument, that returns a named list
  #           with components `loc', `scale' and `shape', each of which is
  #           a numeric m-by-n matrix containing the location, scale and
  #           shape parameters corresponding the the maxima in `data'.
  #
  #     xpar: Covariate information passed to `fpar'.
  #
  #  std.err: Logical. If `TRUE' (the default) then standard errors are
  #           returned as long as the information matrix is invertible.
  #
  #      ...: Optional arguments passed to `optim'.
  #
  # Details:
  #
  #      The block maxima are assumed to be independent and to have GEV
  #      distributions with parameters restricted to belong the family
  #      with form `fpar(p, xpar)' indexed by `p'. The parameterisation is
  #      such that a positive shape parameter corresponds to the Frechet
  #      class.
  #
  # Value:
  #
  #      A list containing the following components:
  #
  #  estimate: A numeric vector with the parameter estimates.
  #
  #   std.err: A numeric vector with the estimated standard errors.
  #
  #       cov: The estimated covariance matrix of the parameter estimates.
  #
  #  deviance: Twice the negative log-likelihood at the maximum-likelihood
  #            estimate.
  #
  # convergence, counts, message: Components taken from the list returned
  #            by 'optim'.
  #
  #       loc: A numeric m-by-n matrix of estimates of the location
  #            parameters corresponding to the maxima in `data'.
  #
  #     scale: A numeric m-by-n matrix of estimates of the scale
  #            parameters corresponding to the maxima in `data'.
  #
  #     shape: A numeric m-by-n matrix of estimates of the shape
  #            parameters corresponding to the maxima in `data'.
  #
  # Warning:
  #
  #      Missing values are ignored, but errors will occur if sufficient
  #      data are missing for parameters to be unidentifiable. Standard
  #      errors assume independence between all block maxima. Users may
  #      use sandwich estimators or bootstrapping to handle dependence.
  #
  # Author:
  #
  #      Chris Ferro <c.a.t.ferro@reading.ac.uk> 14 August 2006
  #
  # References:
  #
  #
  #
  # See Also:
  #
  #
  #
  # Examples:
  #
  #      ns <- 4      # number of sites
  #      nt <- 100    # number of times
  #      library(evd) # data generation
  #
  #      # stationary in space and time
  #
  #      xpar <- c(nt, ns)
  #      fpar <- function(p, xpar) {
  #        loc <- matrix(p[1], xpar[1], xpar[2])
  #        scale <- matrix(p[2], xpar[1], xpar[2])
  #        shape <- matrix(p[3], xpar[1], xpar[2])
  #        list(loc = loc, scale = scale, shape = shape)
  #      }
  #      data <- matrix(rgev(ns * nt, 1, 1, 0), nt, ns)
  #      fit <- fgev.flex(data, c(1, 1, 0), fpar, xpar)
  #      rbind(true = c(1, 1, 0), est = fit$est, ese = fit$std)
  # 
  #      # linear time trend in location parameter; otherwise stationary
  #      
  #      xpar <- matrix(1:nt, nt, ns)
  #      fpar <- function(p, xpar) {
  #        loc <- p[1] + p[2] * xpar
  #        scale <- matrix(p[3], nrow(xpar), ncol(xpar))
  #        shape <- matrix(p[4], nrow(xpar), ncol(xpar))
  #        list(loc = loc, scale = scale, shape = shape)
  #      }
  #      data <- matrix(NA, nt, ns)
  #      for(i in 1:ns) data[, i] <- rgev(nt, 1:nt, 1, 0)
  #      fit <- fgev.flex(data, c(0, 1, 1, 0), fpar, xpar)
  #      rbind(true = c(0, 1, 1, 0), est = fit$est, ese = fit$std)
  # 
  #      # linear spatial trend in log(scale); otherwise stationary
  # 
  #      xpar <- matrix(1:ns, nt, ns, byrow = TRUE)
  #      fpar <- function(p, xpar) {
  #        loc <- matrix(p[1], nrow(xpar), ncol(xpar))
  #        scale <- exp(p[2] + p[3] * xpar)
  #        shape <- matrix(p[4], nrow(xpar), ncol(xpar))
  #        list(loc = loc, scale = scale, shape = shape)
  #      }
  #      data <- matrix(NA, nt, ns)
  #      for(i in 1:ns) data[, i] <- rgev(nt, 0, exp(i), 0)
  #      fit <- fgev.flex(data, c(0, 0, 1, 0), fpar, xpar)
  #      rbind(true = c(0, 0, 1, 0), est = fit$est, ese = fit$std)
  # 
  #      # spatial variation in location; stationary c.v.; else stationary
  # 
  #      xpar <- nt
  #      fpar <- function(p, xpar) {
  #        ns <- length(p) - 2
  #        loc <- matrix(p[1:ns], xpar, ns, byrow = TRUE)
  #        scale <- p[ns+1] * loc
  #        shape <- matrix(p[ns+2], xpar, ns)
  #        list(loc = loc, scale = scale, shape = shape)
  #      }
  #      data <- matrix(NA, nt, ns)
  #      for(i in 1:ns) data[, i] <- rgev(nt, i, 2 * i, 0)
  #      fit <- fgev.flex(data, c(1:ns, 2, 0), fpar, xpar)
  #      rbind(true = c(1:ns, 2, 0), est = fit$est, ese = fit$std)
}

fgev.flex <- function(data, start, fpar, xpar, std.err = TRUE, ...) {
  nll.gev <- function(par) {
    pmat <- fpar(par, xpar)
    loc <- pmat$loc
    scale <- pmat$scale
    shape <- pmat$shape
    if(any(scale <= 0)) return(1e+20)
    gumbel <- (abs(shape) < 1e-06)
    y <- (data - loc) / scale
    z <- 1 + shape * y
    if(any(z <= 0, na.rm = TRUE)) return(1e+20)
    nll <- (1 + 1 / shape) * log(z) + z^(-1 / shape)
    nll[gumbel] <- y[gumbel] + exp(-y[gumbel])
    sum(nll + log(scale), na.rm = TRUE)
  }
  call <- match.call()
  opt <- optim(start, nll.gev, hessian = std.err, ...)
  gev <- fpar(opt$par, xpar)
  out <- list(estimate = opt$par, std.err = rep(NA, length(opt$par)), cov = NULL, deviance = 2 * opt$value, convergence = opt$convergence, counts = opt$counts, message = opt$message, loc = gev$loc, scale = gev$scale, shape = gev$shape)
  cmat <- try(solve(opt$hessian), TRUE)
  
  if(!inherits(cmat, "try-error")) {
    out$std.err <- sqrt(diag(cmat))
    out$cov <- cmat
  }
  structure(c(out, call = call), class = "evd")
}






gevinfom <- function(gamm, scale, n)
  #
  # information matrix for GEV (location, scale, shape) vector
  #
  #
{out = list()
  ceuler <- 0.5772157
  shapeb <- 0.025
  shape <- -gamm
  shape1 <- shape
  if (abs(shape)< 1.e-5) {
    shape <- 1.e-5
  }
  c1 <- 1 - shape
  c2 <- gamma(2 - shape)
  s2 <- scale^2
  sh2 <- shape^2
  p <- c1^2*gamma(1-2*shape)
  q <- c2*(digamma(c1)-c1/shape)
  M11 <- (p)/s2
  M21 <- (p - c2)/(s2*shape)
  M22 <- (1 - 2 * c2 + p)/(s2*sh2)
  M31 <- (q + p/shape)/(scale*shape)
  M32 <- -(1-ceuler-(1-c2)/shape - q - p/shape)/(scale*sh2)
  M33 <- (pi^2/6 + (1-ceuler-1/shape)^2 + 2*q/shape + p/sh2)/sh2
  if (abs(shape)< shapeb) {
    Mneg <- gevinfom(-shapeb*(-1), scale, n)$cov
    Mpos <- gevinfom(shapeb*(-1), scale, n)$cov
    lambda <- (shape+shapeb)/shapeb/2
    M31 <- Mneg[3,1] + (Mpos[3,1]-Mneg[3,1])*lambda
    M32 <- Mneg[3,2] + (Mpos[3,2]-Mneg[3,2])*lambda
    M33 <- Mneg[3,3] + (Mpos[3,3]-Mneg[3,3])*lambda
  }
  M <- n*matrix(c(M11, M21, M31, M21, M22, M32, M31,
                  M32, M33), 3, 3)
  colnames(M) <- c('tail','loc','scale')
  rownames(M) <- c('tail','loc','scale')
  
  cmat <- try(solve(M),TRUE)
  if(!inherits(cmat, "try-error")) {
    out$std.err <- sqrt(diag(cmat))
    out$cov <- cmat 
  }
  
  return(out)
  # index <- match(ci.parameter, c("location", "scale",
  #                                "shape"))
  # sd.param <- sqrt(M.inv[index, index])
}
