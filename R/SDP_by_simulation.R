
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
