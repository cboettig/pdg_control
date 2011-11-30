# file population_models.R
# author Carl Boettiger <cboettig@gmail.com>
# date 2011-11-16
# license BSD
# common population dynamics models with fishing



################################################################
# Library of some possible population dynamics state equations #
################################################################


# temporal variation in f is possible, but increases the memory required, 
# though not the time for the optimization(?)

#' Harvested Beverton Holt growth model
#' @param x fish population 
#' @param h havest level
#' @param p parameters of the growth function, c(A, B), where
#'  A is the maximum growth rate and B the half-maximum in Beverton-Holt.  
#' @returns population next year
#' @details Harvesting takes place before reproduction in ths model.
#'  The carrying capacity is K <- (pars[1]-1)/pars[2]
#'  Try with pars <- c(2,4)
#' @export
BevHolt <- function(x, h, p){
  x <- max(0, x - h)
  A <- p[1] 
  B <- p[2] 
  sapply(x, function(x){ # use sapply so fn accepts vector-valued x
    x <- max(0, x)
    max(0, A * x/(1 + B * x))
  })
}




#' Beverton Holt growth model with fishing-effort based control
#' @param x fish population 
#' @param h fishing effort (harvest is effort times fish population)
#' @param p parameters of the growth function, c(A, B), where
#'  A is the maximum growth rate and B the half-maximum in Beverton-Holt.  
#' @returns population next year
#' @details Harvesting takes place before reproduction in ths model.
#'  The carrying capacity is K <- (pars[1]-1)/pars[2]
#'  Try with pars <- c(2,4)
#' @export
BevHolt_effort <- function(x, h, p){
  S <- max(0, x - x*h) # Escapement from harvest
  A <- p[1] 
  B <- p[2] 
  sapply(S, function(x){ # use sapply so fn accepts vector-valued x
    x <- max(0, x)
    max(0, A * x/(1 + B * x))
  })
}






#' Discrete-time model with an allee effect for alpha > 1
#' @param x the current population level
#' @param h harvest effort
#' @param p vector of parameters c(r, alpha, K) 
#' @returns the population level in the next timestep
#' @details A Beverton-Holt style model with Allee effect.
#'   note that as written, h is fishing EFFORT, not harvest.
#'   Effort above a certain value introduces a fold bifurcation. 
#'   Unharvested carrying capacity is:
#'   K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
#'   The (unharvested) allee theshold is given by:
#'   x = p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 
#'   Bifurcation pt is h = (p[1]*sqrt(p[3])-2)/2 
#'   Try with pars = c(1,2,6), h=.01
#' 
#'   Consider updating to be a function of x-h, instead?
#' @export
Myer <- function(x, h, p){
   max(0, p[1] * x ^ p[2] / (1 + x ^ p[2] / p[3])  - h * x)
}

#' Discrete-time model with an allee effect for alpha > 1 with harvest control
#' @param x the current population level
#' @param h total harvest level
#' @param p vector of parameters c(r, alpha, K) 
#' @returns the population level in the next timestep
#' @details A Beverton-Holt style model with Allee effect.
#'   Unharvested carrying capacity is:
#'   K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
#'   The (unharvested) allee theshold is given by:
#'   x = p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 
#' @export
Myer_harvest <- function(x, h, p){
   max(0, p[1] * x ^ p[2] / (1 + x ^ p[2] / p[3])  - h)
}






#' Ricker-like model with Allee effect (Allen)
#' @param x the current population level
#' @param h harvest level 
#' @param p a vector of parameters c(r, K, C) 
#' @returns the population level in the next timestep
#' @export
RickerAllee <- function(x, h, p){
    x <- max(0,x-h)
    x * exp(p[1] * (1 - x / p[2]) * (x - p[3]) / p[2] ) 
}


#' Basic Ricker model 
#' @param x the current population level
#' @param h harvest level 
#' @param p a vector of parameters c(r, K) 
#' @returns the population level in the next timestep
#' @export
Ricker <- function(x,h,p){
  x <- max(0, x-h) 
  max(0, x * exp(p[1] * (1 - x / p[2] )) )
}


#' Coral-Parrotfish model
#' @param x vector of population levels: macroalgae, coral, parrotfish
#' @param h harvesting effort on parrotfish
#' @param p c(a, g, T, gamma, r, d, s, K)
#'            1  2  3     4   5  6  7  8
#' @details 
#' @export
coral <- function(x, h, p){
 x_t1 <- p[1] * x[1] * x[2] + p[2] * x[3] * x[1] / (x[1] + p[3]) + p[4] * p[3] * x[1]
 x_t2 <- (p[5] * p[3] - p[6] - p[1] * x[1]  ) * x[2] 
 x_t3 <- p[7] * x[3] * (1 - x[3] / (p[8] * x[2])) - h * x[3]
 c(x_t1, x_t2, x_t3)
}







