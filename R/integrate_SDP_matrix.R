
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