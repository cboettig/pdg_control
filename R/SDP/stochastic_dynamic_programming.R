# file stochastic_dynamic_programming.R
# author Carl Boettiger <cboettig@gmail.com>
# date 2011-11-07
# licence BSD
# adapted from SDP.m by Michael Bode
# 
# Contains functions for stochastic dynamic programming 
# see Reed.R for example uses 

#########################################################################
# A function to generate the transition matrix used for the SDP routine #
#########################################################################

#' Determine the Stochastic Dynamic Programming matrix.
#' @param f the growth function of the escapement population (x-h)
#'   should be a function of f(t, y, p), with parameters p
#' @param p the parameters of the growth function
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param sigma the variance of the noise process
#' @returns the transition matrix at each value of h in the grid.  
determine_SDP_matrix <- function(f, p, x_grid, h_grid, sigma){
  SDP_Mat <- lapply(h_grid, function(h){
    SDP_matrix <- matrix(0, nrow=gridsize, ncol=gridsize)
    # Cycle over x values
    for(i in 1:gridsize){ ## VECTORIZE ME
      ## Calculate the 
      x1 <- x_grid[i]
      x2_expected <- f(x1, h, pars)
      ## If expected 0, go to 0 with probabilty 1
      if( x2_expected == 0) 
        SDP_matrix[i,1] <- 1  
      else {
        # relative probability of a transition to that state
        ProportionalChance <- x_grid / x2_expected
        # lognormal due to multiplicative Gaussian noise
        Prob <- dlnorm(ProportionalChance, 0, sigma)
        # Store normalized probabilities in row
        SDP_matrix[i,] <- Prob/sum(Prob)
      }
    }
    SDP_matrix
  })
  SDP_Mat
}


########################################################################
# A function to identify the dynamic optimum using backward iteration  #
########################################################################

#' Identify the dynamic optimum using backward iteration (dynamic programming)
#' @param SDP_Mat the stochastic transition matrix at each h value
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param OptTime the stopping time 
#' @param xT the boundary condition population size at OptTime
#' @param c the cost/profit function, a function of harvested level
#' @param delta the exponential discounting rate
#' @returns list containing the matrices D and V.  D is an x_grid by OptTime
#'  matrix with the indices of h_grid giving the optimal h at each value x as 
#'  the columns, with a column for each time.  
#'  V is a matrix of x_grid by x_grid, which is used to store the value 
#'  function at each point along the grid at each point in time.  
#'  The returned V gives the value matrix at the first (last) time.  
find_dp_optim <- function(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, 
                          delta, reward=10){

  gridsize <- length(x_grid)
  HL <- length(h_grid)

  ## Initialize space for the matrices
  D <- matrix(NA, nrow=gridsize, ncol=OptTime)
  V <- rep(0,gridsize) # initialize BC, 
  # give a fixed reward for having value larger than xT at the end. 
  V[ x_grid >= xT ] <- reward # "Scrap Value" for x(T) >= xT
  # loop through time  
  for(time in 1:OptTime){
    # try all potential havest rates
    V1 <- sapply(1:HL, function(i){
      # Transition matrix times V gives dist in next time
      SDP_Mat[[i]] %*% V + 
      # then (add) harvested amount times discount
       profit(x_grid, h_grid[i]) * exp(-delta * (OptTime - time))
    })

    # find havest, h that gives the maximum value
    out <- sapply(1:gridsize, function(j){
      value <- max(V1[j,], na.rm = T) # each column is a diff h, max over these
      index <- which.max(V1[j,])      # record index so we can recover the h's 
      c(value, index)                 # returns both values 
    })

    # V{t+1} = max_h V{t} at each possible state value, x
    V <- out[1,]                        # The new value-to-go
    D[,OptTime-time+1] <- out[2,]       # The index positions
  }

  # Reed derives a const escapement policy saying to fish the pop down to:
  ReedThreshold <- x_grid[sum(D[,1]==1)] # easy way
  # calculation is harder for general f, need to start at top
  # finds the largest population for which you shouldn't harvest: 
  ReedThreshold <- x_grid[max(which(D[,1] == 1))]
  list(D=D, V=V, S=ReedThreshold)
}



#' Forward simulate given the optimal havesting policy, D
#' @param f the growth function of the escapement population (x-h)
#'   should be a function of f(y, p), with parameters p
#' @param p the parameters of the growth function
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param sigma the variance of the noise process
#' @param Xo initial stock size
#' @param D the optimal solution indices on h_grid, 
#'  given for each possible state at each timestep
#' @param sigma_assess amount of uncertainty in the assessment of stock size
#' @param sigma_harvest amount of noise in implementing quotas 
#' @param interval is the years between updating the harvest quota
#' @returns a data frame with the time, fishstock, harvested amount,
#'  and what the stock would have been without that year's harvest.  
ForwardSimulate <- function(f, pars, x_grid, h_grid, sigma, x0, D,
                            sigma_assess = 0, sigma_harvest = 0,
                            interval = 1){
  # shorthand names
  n <- x_grid
  h <- h_grid

  # initialize variables with initial conditions
  OptTime <- dim(D)[2]
  gridsize <- length(n)
  x_h <- numeric(OptTime) # population dynamics with harvest
  x   <- numeric(OptTime) # What would happen with no havest
  h   <- numeric(OptTime) # optimal havest level
  x_h[1] <- x0 
  x[1]   <- x0 

  for(t in 1:(OptTime-1)){
    # Current state; can add noise to stock assessment, x_h[t] 
    St <- which.min(abs(n - x_h[t] * rnorm(1, 1, sigma_assess))) 
    # Set harvest quota on update years
  # could update harvest quota annually based on projections, but assess stock only periodically.  Could update stock annually, but adjust quota periodically.  
    if(t %% interval == 0){ 
      h[t + 1] <- h_grid[D[St,t + 1]] 
    } else {
      h[t + 1] <- h[t]
    }
    h[t + 1] <- h[t + 1] * rnorm(1, 1, sigma_harvest) 
    z <- rnorm(1,1,sigma)            # Noise in growth
    x_h[t+1] <- z*f(x_h[t], h[t + 1], pars) # with havest
    x[t+1]   <- z * f(x[t], 0, pars) # havest-free dynamics
  }
  data.frame(time=1:OptTime, fishstock=x_h, harvest=h, unharvested=x) 
}


################################################################
# Library of some possible population dynamics state equations #
################################################################


# temporal variation in f is possible, but increases the memory required, 
# though not the time for the optimization(?)

####### Define our population dynamics / state equation  #######
#' Harvested Beverton Holt growth model
#' @param x fish population 
#' @param h havest level
#' @param p parameters of the growth function, c(A, B), where
#'  A is the maximum growth rate and B the half-maximum in Beverton-Holt.  
#' @returns population next year
#' @details Harvesting takes place before reproduction in ths model.
#'  The carrying capacity is K <- (pars[1]-1)/pars[2]
#'  Try with pars <- c(2,4)
BevHolt <- function(x, h, p){
  x <- max(0, x - h)
  A <- p[1] 
  B <- p[2] 
  sapply(x, function(x){ # use sapply so fn accepts vector-valued x
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
Myer <- function(x, h, p){
   max(0, p[1] * x ^ p[2] / (1 + x ^ p[2] / p[3])  - h * x)
}

#' Ricker-like model with Allee effect (Allen)
#' @param x the current population level
#' @param h harvest level 
#' @param p a vector of parameters c(r, K, C) 
#' @returns the population level in the next timestep
RickerAllee <- function(x, h, p){
    x <- max(0,x-h)
    x * exp(p[1] * (1 - x / p[2]) * (x - p[3]) / p[2] ) 
}


Ricker <- function(x,h,p){
  x <- max(0, x-h) 
  max(0, x * exp(p[1] * (1 - x / p[2] )) )
}

#' corals

#' Coral-Parrotfish model
#' @param x vector of population levels: macroalgae, coral, parrotfish
#' @param h harvesting effort on parrotfish
#' @param p c(a, g, T, gamma, r, d, s, K)
#'            1  2  3     4   5  6  7  8
#' @details 
coral <- function(x, h, p){
 x_t1 <- p[1] * x[1] * x[2] + p[2] * x[3] * x[1] / (x[1] + p[3]) + p[4] * p[3] * x[1]
 x_t2 <- (p[5] * p[3] - p[6] - p[1] * x[1]  ) * x[2] 
 x_t3 <- p[7] * x[3] * (1 - x[3] / (p[8] * x[2])) - h * x[3]
 c(x_t1, x_t2, x_t3)
}






