# file: training_prob2_collocation.R
# author: Carl Boettiger, \url{http://carlboettiger.info}
# date: 2011-11-24

# Libraries
require(bvpSolve)

#' State Equation(s): Fish population dynamics
#' @param t is time 
#' @param y is a vector of (x,h)', the fish pop and harvest level
#' @param p parameters, c(alpha, K, C) 
#' @return \dot x = f(x), population growth
f <- function(t, y, p){
  # Rename explicitly so equation is easier to read but still fast.
  x <- y[1]
  h <- y[2]
  alpha <- pars[1]
  K <- pars[2]
  C <- pars[3]
  x * alpha * ((K - x) / K) * ((x - C) / K) - h * x
}
#' Derivative with respect to state x
df <- function(t, y, p){
  x <- y[1]
  h <- y[2]
  alpha <- pars[1]
  K <- pars[2]
  C <- pars[3]
  - alpha * ( C* K - 2 * x * K - 2 * x * C + 3 * x ^ 2) / K ^ 2 - h
}

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


