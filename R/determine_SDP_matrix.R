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




