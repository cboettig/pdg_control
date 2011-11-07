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
# Fish population dynamics:
# X_{t+1} = Z_n f(X_n) 
# f(x) = 
# 
# 


require(ggplot2) # nicer plotting package
rm(list=ls())   # Start wtih clean workspace 


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
      x2_expected <- f(x1-h, pars)
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
find_dp_optim <- function(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta){
  gridsize <- length(x_grid)
  HL <- length(h_grid)

  ## Initialize space for the matrices
  D <- matrix(NA, nrow=gridsize, ncol=OptTime)
  V <- rep(0,gridsize) # initialize BC, "No scrap value" (leave no fish)

# loop through time  
  for(time in 1:OptTime){
    # try all potential havest rates
    V1 <- sapply(1:HL, function(i){

      # havest cannot exceed population size
      min_hn <- sapply(x_grid, function(n) min(h_grid[i], n))

      # Transition matrix times V gives dist in next time
      # then (add) harvested amount times discount
      SDP_Mat[[i]] %*% V + profit(min_hn) * exp(-delta * (OptTime-time))
    })

    # find havest, h that gives the maximum value
    out <- sapply(1:gridsize, function(j){
      value <- max(V1[j,], na.rm = T)   # each column is a diff h, max over these
      index <- which.max(V1[j,])        # record index so we can recover the h's 
      c(value, index)                   # returns both values 
    })

    # V{t+1} = max_h V{t} at each possible state value, x
    V <- out[1,]                        # The new value-to-go
    D[,OptTime-time+1] <- out[2,]       # The index positions
  }
  list(D=D, V=V)
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
#' @returns a data frame with the time, fishstock, harvested amount,
#'  and what the stock would have been without that year's harvest.  
ForwardSimulate <- function(f, pars, x_grid, h_grid, sigma, x0, D){
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
    St <- which.min(abs(n - x_h[t])) # Current state
    h[t] <- h_grid[D[St,t+1]]      # Optimal harvest for state
    z <- rnorm(1,1,sigma)
    x_h[t+1] <- z*f(x_h[t]-h[t], pars) # with havest
    x[t+1]   <- z*f(x_h[t], pars) # no havest
  }
  data.frame(time=1:OptTime, fishstock=x_h, harvest=h, unharvested=x) 
}




#############################################################################
#                                                                           #
#  These are probably the functions you want to modify to specify a         #
#  particular problem.  The above functions are somewhat general and        #
#  used to implement the SDP solution on the system defined below.          #
#                                                                           #
#############################################################################



########## Define our population dynamics / state equation  ################ 
#' Harvested Beverton Holt growth model
#' @param x fish population that reproduces (usually x-h)
#' @param p parameters of the growth function, c(A, B), where
#'  A is the maximum growth rate and B the half-maximum in Beverton-Holt.  
#' @returns population next year
BevHolt <- function(x, p){
  A <- p[1] 
  B <- p[2] 
  sapply(x, function(x){ # use sapply so fn accepts vector-valued x
    x <- max(0, x)
    max(0, A*x/(1+B*x))
  })
}



# Define all parameters 
delta <- 0.1      # economic discounting rate
OptTime <- 25     # stopping time
sigma <- 0.2      # Noise process
gridsize <- 5     # gridsize (discretized population)
pars <- c(2,4)    # Beverton-Holt/f(x) pars, A, B
K <- (pars[1]-1)/pars[2]   # Unharvested deterministic equib pop


# define a profit function, price minus cost
profit <- function(h){
  p <- 1
  c <- 0.001 # higher extraction costs result in less fishing 
  # sapply: support for vector-valued h
  sapply(h, function(h) max(0,p*h - c/h))
}

# Set up the grid 
x_grid <- seq(0, 2*K, length=gridsize)  # population size
h_grid <- x_grid  # vector of havest levels, use same res as stock

# Calculate the transition matrix 
SDP_Mat <- determine_SDP_matrix(BevHolt, pars, x_grid, h_grid, sigma)

# Find the optimum by dynamic programming 
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, 0, profit, delta)

# Simulate the optimal routine on a stoch realization of growth dynamics
out <- ForwardSimulate(BevHolt, pars, x_grid, h_grid, sigma, K/2, opt$D)

# Plot the results 
ggplot(out) + geom_line(aes(time, fishstock)) + 
  geom_line(aes(time, unharvested), col=I("green"))









## as area plot
# data <- melt(data, id="time")
# ggplot(data, aes(year, value, fill=variable)) + geom_area()


#ReedThreshold <- n[ sum(D[,1]==1) ]

