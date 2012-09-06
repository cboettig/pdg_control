# file determine_SDP_matrix.R 
# author Carl Boettiger <cboettig@gmail.com>
# date 2011-11-07
# licence BSD
# SDP solution adapted from SDP.m by Michael Bode
# 
# Contains three different functions for determining the 
# Stochastic Transition Matrix used in the Dynamic Programming solution
#


#########################################################################
# A function to generate the transition matrix used for the SDP routine #
#########################################################################

#' Determine the Stochastic Dynamic Programming matrix.
#' @param f the growth function of the escapement population (x-h)
#'   should be a function of f(t, y, p), with parameters p
#' @param p the parameters of the growth function
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param sigma_g the variance of the population growth process
#' @param pdfn the probability density function, taking the proportional
#' chance of a transition to that state, f(x), and the parameter sigma_g
#' By default it will use the log normal density.  
#' @return the transition matrix at each value of h in the grid. 
#' @details this analytical approach doesn't reliably support other 
#'  sources of variation.  The quality of the analytic approximations 
#'  (lognormal) can be tested. 
#' @export
determine_SDP_matrix <- function(f, p, x_grid, h_grid, sigma_g,
                                 pdfn=function(P, s) dlnorm(P, 0, s)){
  gridsize <- length(x_grid)
  SDP_Mat <- lapply(h_grid, function(h){
    SDP_matrix <- matrix(0, nrow=gridsize, ncol=gridsize)
    # Cycle over x values
    for(i in 1:gridsize){ ## VECTORIZE ME
      ## Calculate the expected transition  
      x1 <- x_grid[i]
      x2_expected <- f(x1, h, p)
      ## If expected 0, go to 0 with probabilty 1
      if( x2_expected == 0) 
        SDP_matrix[i,1] <- 1  
      else {
        # relative probability of a transition to that state
        ProportionalChance <- x_grid / x2_expected
        Prob <- pdfn(ProportionalChance, sigma_g)
        # Store normalized probabilities in row
        SDP_matrix[i,] <- Prob/sum(Prob)
      }
    }
    SDP_matrix
  })
  SDP_Mat
}


#' Determine the Stochastic Dynamic Programming matrix.
#'
#' Integrate the multidimensional function using the cubature package
#' @param f the growth function of the escapement population (x-h)
#'   should be a function of f(t, y, p), with parameters p
#' @param p the parameters of the growth function
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param sigma_g the variance of the population growth process
#' @param sigma_m noise in stock assessment (currently assumes lognormal)
#' @param sigma_i noise in implementing the quota (lognormal)
#' @return the transition matrix at each value of h in the grid. 
#' @import cubature
#' @export
integrate_SDP_matrix  <- function(f, p, x_grid, h_grid, sigma_g, sigma_m, sigma_i){
  gridsize <- length(x_grid)
  SDP_Mat <- lapply(h_grid, function(h){
    mat <- sapply(x_grid, function(y){
      # Handle the case of 0 expectation seperately, maps 100% to extinction
      expected <- f(y,h,p)
      if(expected==0){
        Prob <- numeric(gridsize)
        Prob[1] <- 1
      } else {
        # dividing x by the expected value is same as scaling distribution to mean 1
        pdf_zg <- function(x, expected) dlnorm(x/expected, 0, sigma_g)
        pdf_zm <- function(x) dlnorm(x, 0, sigma_m)
        pdf_zi <- function(x,q) dlnorm(x, log(q), sigma_i)
        Prob <- sapply(x_grid, function(y){
          F <- function(x) 
            pdf_zg(y, f(x[1], x[2], p)) * pdf_zm(x[1]) * pdf_zi(x[2], h)
          int <- adaptIntegrate(F, c(0, 0), c(10*K, 10*K))
          int$integral
        })   
      }
      Prob/sum(Prob)
    })
  t(mat)
  })
  SDP_Mat
}


#' Determine the transtion matrix using stochastic simulation
#' @param f the growth function of the escapement population (x-h)
#'   should be a function of f(t, y, p), with parameters p
#' @param p the parameters of the growth function
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param z_g a function determing the probability distrubtion for the 
#'  stochastic population growth process (draws a random variable z_g)
#' @param z_m a function determining the porbability distribution for
#'  measurement uncertainty in the assessment of stock size (random variable)
#' @param z_i function for implementation uncertainty in quotas 
#' @return the transition matrix at each value of h in the grid.  
#' @import snowfall
#' @import ggplot2
#' @import Hmisc
#' @export
SDP_by_simulation <- function(f, p, x_grid, h_grid, z_g, z_m, z_i, reps = 999){
  require(snowfall) # support parallelization of this
  sfExportAll()
  sfLibrary(ggplot2) # for the bin function 

  bw <- x_grid[2] - x_grid[1]
  lower <-x_grid[1]
  upper <- x_grid[length(x_grid)] + bw

  SDP_Mat <- sfLapply(h_grid, function(h){ 
    mat <- sapply(x_grid, function(x){
      x_t1 <- replicate(reps, z_g() * z_m() * f(x / z_m(), z_i() * h, p)) 
      a <- mybin( x_t1, binwidth=bw, range=c(lower, upper))$count
      a / sum(a) 
    })
    t(mat)
  })
  SDP_Mat
}





# Bin data from old version of ggplot2
# This function powers \code{\link{stat_bin}}R
#
# @keyword internal
mybin <- function(x, weight=NULL, binwidth=NULL, origin=NULL, breaks=NULL, range=NULL, width=0.9, drop = FALSE) {
  
  if (is.null(weight))  weight <- rep(1, length(x))
  weight[is.na(weight)] <- 0

  if (is.null(range))    range <- range(x, na.rm = TRUE, finite=TRUE)
  if (is.null(binwidth)) binwidth <- diff(range) / 30

  if (is.integer(x)) {
    bins <- x
    x <- sort(unique(bins))
    width <- width    
  } else if (diff(range) == 0) {
    width <- width
    bins <- x
  } else { # if (is.numeric(x)) 
    if (is.null(breaks)) {
      if (is.null(origin)) {
        breaks <- fullseq(range, binwidth)        
      } else {
        breaks <- seq(origin, max(range) + binwidth, binwidth)
      }
    }
    bins <- cut(x, sort(breaks), include.lowest=TRUE)
    left <- breaks[-length(breaks)]
    right <- breaks[-1]
    x <- (left + right)/2
    width <- diff(breaks)
  }

  results <- data.frame(
    count = as.numeric(tapply(weight, bins, sum, na.rm=TRUE)),
    x = x,
    width = width
  )

  res <- within(results, {
    count[is.na(count)] <- 0
    density <- count / width / sum(count, na.rm=TRUE)
    ncount <- count / max(count, na.rm=TRUE)
    ndensity <- density / max(density, na.rm=TRUE)
  })
  if (drop) res <- subset(res, count > 0)
  res
}

# Generate sequence of fixed size intervals covering range
# All locations are multiples of size
# 
# @arguments range
# @arguments interval size
# @keyword internal
# @seealso \code{\link{reshape}{round_any}}
fullseq <- function(range, size) {
  seq(
    plyr::round_any(range[1], size, floor), 
    plyr::round_any(range[2], size, ceiling), 
    by=size
  )
}

