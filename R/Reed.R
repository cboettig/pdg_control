# file Reed.R
# author Carl Boettiger, <cboettig@gmail.com>
# date 2011-11-02
# modified from SDP.m, by Michael Bode.  
# 
# Implements a numerical version of the SDP described in:
#   Reed, W.J., 1979. Optimal Escapement Levels in Stochastic
#   and Deterministic Harvesting Models. Journal of Environmental 
#   Economics and Management. 6: 350-363.
#
#
# 
# 
# Fish population dynamics:
# X_{t+1} = Z_n f(X_n) 
# f(x) = 



rm(list=ls())   # Start wtih clean workspace 

############## Parameters ################
delta <- 0.1      # economic discounting rate
sigma <- 0.4      # Noise process
gridsize <- 100   # gridsize (discretized population)
A <- 2            # Beverton-Holt/f(x) pars
B <- 4            # Beverton-Holt/f(x) pars
K <- (A-1)/B      # Unharvested deterministic equilib population



#' Harvested Beverton Holt growth model
#' @param x fish population currently
#' @param A growth rate 
#' @param B half-maximum
#' @param h havest
#' @returns population next year
SubBevHolt <- function(x1, A, B, h){
  sapply(x1, function(x){ # accept vector-valued x
    x_minus_h <- max(0, x-h)
    max(0, A*x_minus_h/(1+B*x_minus_h))
  })
}

# Show the population growth function f(x)
#curve(SubBevHolt(x,A,B,0), 0, 2*K)

n_vec <- seq(0, 2*K, length=gridsize)  # population size
d <- n_vec[2] - n_vec[1]               # grid spacing
H_vec <- n_vec                         # vector of havest levels

# Set up the grid
n_vec <- seq(0, 2*K, length=gridsize)  # population size
H_vec <- n_vec                         # vector of havest levels

#' Caculate the Transition matrix 
#' @details For each harvest level, generate a transition
#' matrix for the probability of going form any population 
#' level to any other population level in the grid at time t.
SDP_Mat <- lapply(H_vec, function(h){
  SDP_matrix <- matrix(0, nrow=gridsize, ncol=gridsize)
  # Cycle over x values
  for(i in 1:gridsize){
    ## Calculate the 
    x1 <- n_vec[i]
    x2_expected <- SubBevHolt(x1, A, B, h)
    ## If expected 0, go to 0 with probabilty 1
    if( x2_expected == 0) 
      SDP_matrix[i,1] <- 1  
    else {
      ProportionalChance <- n_vec / x2_expected
      Prob <- dlnorm(ProportionalChance, 0, sigma)
      SDP_matrix[i,] <- Prob/sum(Prob)
    }
  }
  SDP_matrix
})

##########################################################
# Identify the dynamic optimum using backward iteration  #
##########################################################
OptTime = 25
D <- matrix(NA, nrow=gridsize, ncol=OptTime)
V <- rep(0,gridsize) # initialize, "No scrap value" (leave no fish)
V <- rep(.1,gridsize) # initialize, 

# loop through time  
for(time in 1:OptTime){

  # try all potential havest rates
  V1 <- sapply(1:gridsize, function(i){

    # havest cannot exceed population size
    min_hn <- sapply(n_vec, function(n) min(H_vec[i], n))

    # Transition matrix times V gives dist in next time
    # then (add) harvested amount times discount
    SDP_Mat[[i]] %*% V + min_hn * exp(-delta * (OptTime-time))
  })

  # find havest, h that gives the maximum value 
  out <- sapply(1:gridsize, function(j){
    value <- max(V1[j,], na.rm = T) 
    index <- which.max(V1[j,])
    c(value, index)
  })

  V <- out[1,]                  # The new value-to-go
  D[,OptTime-time+1] <- out[2,] # The index positions
}




############
#  Plots 
############






#' Forward simulate given the optimal havesting policy, D
#' @param Xo inital stock size 
#' @param 
ForwardSimulate <- function(x0, A, B, D, sigma, n, H){

# initialize variables with initial conditions
  OptTime <- dim(D)[2]
  gridsize <- length(n)
  x_h <- numeric(OptTime) # population dynamics with harvest
  x   <- numeric(OptTime) # What would happen with no havest
  h   <- numeric(OptTime) # optimal havest level
  x_h[1] <- x0 
  x[1]   <- x0 

  for(t in 1:(OptTime-1)){
    St <- which.min(abs(n - x_h[t])) # Current state
    h[t] <- H[D[St,t+1]]      # Optimal harvest for state
    z <- rnorm(1,1,sigma)
    x_h[t+1] <- z*SubBevHolt(x_h[t], A, B, h[t]) # with havest
    x[t+1]   <- z*SubBevHolt(x_h[t], A, B, 0) # no havest
  }
  ReedThreshold <- n[ sum(D[,1]==1) ]
  data.frame(time=1:OptTime, fishstock=x_h, harvest=h, unharvested=x) 
}

out <- ForwardSimulate(K/2, A, B, D, sigma, n_vec, H_vec)

require(ggplot2)
qplot(time, fishstock, data=out)

