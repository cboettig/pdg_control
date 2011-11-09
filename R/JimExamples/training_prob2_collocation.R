# file: training_prob2_collocation.R
# author: Carl Boettiger, \url{http://carlboettiger.info}
# date: 2011-11-24

# Libraries
require(bvpSolve)

##################################################################
# Solve the ODE system of fish dynamics at fixed harvest level   #
##################################################################
fish <- function(t,y,p){
  dy1 <- f(t, y, p) 
  dy2 <- 0                # constant harvest level 
  list(c(dy1, dy2), NULL)
}
jac <- function(t,y,p){
  matrix(c(df(t, y, p), 0,
           0,           0),
         2,2, byrow=T)
}

##################################################################
# Specify & solve the Boundary Value Problem                     #
##################################################################

#' fun defines the bvp system to be solved
#' @param t time variable 
#' @param y a vector of the system state y[1],y[2] = (x,h)
#' @param p parameters: c(alpha, K, C, gamma, rho)
fun <- function(t,y,p){
  gamma <- p[4]
  rho <- p[5]

  dy1 <- f(t, y, p) 
  dy2 <- (1 / gamma) * (rho - df(t,y,p))

  list(c(dy1, dy2))
}


