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
#' @return population next year
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
#' @return population next year
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
#' @return the population level in the next timestep
#' @details Unharvested carrying capacity is:
#'   K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
#'   The (unharvested) allee theshold is given by:
#'   x = p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 
#'   Bifurcation pt is h = (p[1]*sqrt(p[3])-2)/2 
#'   Try with pars = c(1,2,6), h=.01
#' 
#'   Consider updating to be a function of x-h, instead?
#' @export
Myers <- function(x, h, p){
 sapply(x, function(x){
   x <- max(0, x - h) 
  max(0, p[1] * x ^ p[2] / (1 + x ^ p[2] / p[3]) )
 })
}

#' Discrete-time model with an allee effect for alpha > 1 with harvest control
#' @param x the current population level
#' @param h total harvest level
#' @param p vector of parameters c(r, alpha, K) 
#' @return the population level in the next timestep
#' @details A Beverton-Holt style model with Allee effect.
#'   Unharvested carrying capacity is:
#'   K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
#'   The (unharvested) allee theshold is given by:
#'   x = p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 
#' @export
Myer_harvest <- function(x, h, p){
   sapply(x, function(x) max(0, p[1] * x ^ p[2] / (1 + x ^ p[2] / p[3])  - h))
}




#' Harvesting model with a bifurcation from May (1977)
#' increasing parameter "a" (p[3]) will drive the bifurcation
#' @param x the current population level
#' @param h total harvest level
#' @param p vector of parameters c(r, K, a, H, Q) 
#' @return the population level in the next timestep
#' @details try pars <- c(r = .75, k = 10, a=1, H=1, Q = 3)
#' See the position of the bifurcation at around a = 1.9:
#' curve(.75*(1-x/10), 0, 10)
#' curve(1.9*x^2/(x^3+1), 0, 10, add=T, col="red")
#' goes to nonzero alternate stable state.  
#' @export
May <- function(x, h, p){
  sapply(x, function(x){
         s <- max(x - h, 0) # escapement
         r <- as.numeric(p[1])
         K <- as.numeric(p[2])
         a <- as.numeric(p[3])
         H <- as.numeric(p[4])
         Q <- as.numeric(p[5])
         s * exp(r * (1 - s / K) - a * s ^ (Q - 1) / (s ^ Q + H ^ Q)) 
  })
}


#' Ricker-like model with Allee effect (Allen)
#' @param x the current population level
#' @param h harvest level 
#' @param p a vector of parameters c(r, K, C) 
#' @return the population level in the next timestep
#' @export
RickerAllee <- function(x, h, p){
  sapply(x, function(x){ 
    x <- max(0,x-h)
    x * exp(p[1] * (1 - x / p[2]) * (x - p[3]) / p[2] ) 
  })
}


#' Basic Ricker model 
#' @param x the current population level
#' @param h harvest level 
#' @param p a vector of parameters c(r, K) 
#' @return the population level in the next timestep
#' @export
Ricker <- function(x,h,p){
  sapply(x, function(x){ 
    x <- max(0, x-h) 
    max(0, x * exp(p[1] * (1 - x / p[2] )) )
  })
}


#' Coral-Parrotfish model
#' @param x vector of population levels: macroalgae, coral, parrotfish
#' @param h harvesting effort on parrotfish
#' @param p c(a, g, T, gamma, r, d, s, K)
#'            1  2  3     4   5  6  7  8
#' @references Blackwood et al. (2011) doi:10.1890/10-2195.1
#' @details 
#' $$\begin{align}
#' \frac{dM}{dt} = aMC - \frac{g(P) M}{M+T} + \gamma M T \\
#' \frac{dC}{dt} = rTC - dC - a M C \\
#' \frac{dP}{dt} = sP \left( 1- \frac{P}{\beta K(C) } \right) - h P 
#' \end{align}$$
#' @export
coral <- function(x, h, p = c(a = 0.1, g = 1, T = , gamma = 0.8, r = 1, d = 0.44, s= 0.49, K = 1)){
 x_t1 <- p[1] * x[1] * x[2] + p[2] * x[3] * x[1] / (x[1] + p[3]) + p[4] * p[3] * x[1]
 x_t2 <- (p[5] * p[3] - p[6] - p[1] * x[1]  ) * x[2] 
 x_t3 <- p[7] * x[3] * (1 - x[3] / (p[8] * x[2])) - h * x[3]
 c(x_t1, x_t2, x_t3)
}
                  


#' Coral-Parrotfish model
#' @param x vector of population levels: macroalgae, coral, parrotfish
#' @param h harvesting effort on parrotfish
#' @param p c(a, g, gamma, r, d)
#' @param dt descrization scale
#' @references Mumby et al. (2007) doi:10.1038/nature06252
#' Blackwood et al. (2011) doi:10.1890/10-2195.1
#' @details 
#' $$\begin{align}
#' \frac{dM}{dt} = aMC - \frac{(g-h) M}{M+T} + \gamma M T \\
#' \frac{dC}{dt} = rTC - dC - a M C \\
#' T = 1 - M - C
#' \end{align}$$
#' @export
mumby <- function(x, h, p= c(a = 0.1, g = .5, gamma = 0.8, r = 1, d = 0.44), dt=0.025){
  M <- x[1] # Macroalgae
  C <- x[2] # Corals
  a <- p[1] # algal growth rate
  g <- p[2] # grazing rate (reduced by harvesting)
  gamma <- p[3] # algal colinization of dead tufts
  r <- p[4] # Coral growth rate
  d <- p[5] # Coral death rate
  M_t <- M * exp( dt * (a * M * C - (g - h) * M / (M + (1-M-C)) + gamma * M * (1-M-C)) )
  C_t <- C * exp( dt * (r * (1-M-C) * C - d * C - a * M * C) )
  c(M_t, C_t)
}





