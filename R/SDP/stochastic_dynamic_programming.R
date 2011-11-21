# file stochastic_dynamic_programming.R
# author Carl Boettiger <cboettig@gmail.com>
# date 2011-11-07
# licence BSD
# adapted from SDP.m by Michael Bode
# 
# Contains functions for stochastic dynamic programming 
# see Reed.R for example uses 




#########################################################################
# A library of distribution functions for sources of stohcasticity      #
#########################################################################

# z_g is the stochasticity in the growth process, x_{t+1} = z_g f(x_t)
# z_m is the measurement error in the stock assessment, m_t = z_m x_t
# z_i is the implementation error in the harvest quota: h_t = z_i q_t

# Normal random vars -- an unusual choice given the negative domain support
z_g <- function(sigma_g) rnorm(1,1, sigma_g)
z_m <- function(sigma_m) rnorm(1,1, sigma_m)
z_i <- function(sigma_i) rnorm(1,1, sigma_i)

# Log-normal distribution -- perhaps the most natural, at least for z_g
# mean is 1 = exp(mu + sigma^2/2), then
# exp(1) - sigma^2/2 = mu
z_g <- function() rlnorm(1,  log(1)-sigma_g^2/2, sigma_g) # mean 1




# No noise
#z_g <- function() 1
z_m <- function() 1
z_i <- function() 1


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
#' @returns the transition matrix at each value of h in the grid. 
#' @details this analytical approach doesn't reliably support other 
#'  sources of variation.  The quality of the analytic approximations 
#'  (lognormal) can be tested.  
determine_SDP_matrix <- function(f, p, x_grid, h_grid, sigma_g){
  gridsize <- length(x_grid)
  SDP_Mat <- lapply(h_grid, function(h){
    SDP_matrix <- matrix(0, nrow=gridsize, ncol=gridsize)
    # Cycle over x values
    for(i in 1:gridsize){ ## VECTORIZE ME
      ## Calculate the 
      x1 <- x_grid[i]
      x2_expected <- f(x1, h, p)
      ## If expected 0, go to 0 with probabilty 1
      if( x2_expected == 0) 
        SDP_matrix[i,1] <- 1  
      else {
        # relative probability of a transition to that state
        ProportionalChance <- x_grid / x2_expected
        # lognormal due to multiplicative Gaussian noise
        Prob <- dlnorm(ProportionalChance, 0, sigma_g)
        # Store normalized probabilities in row
        SDP_matrix[i,] <- Prob/sum(Prob)


## dnorm(x_i, mu, sigma)  x_i+1,  
      }
    }
    SDP_matrix
  })
  SDP_Mat
}


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
#' @returns the transition matrix at each value of h in the grid.  
#' @import snowfall
#' @import ggplot2
#' @export
SDP_by_simulation <- function(f, p, x_grid, h_grid, z_g, z_m, z_i, reps = 999){
  require(snowfall) # support parallelization of this
  sfExportAll()
  sfLibrary(ggplot2) # for the bin function 
  SDP_Mat <- sfLapply(h_grid, function(h){ 
    mat <- sapply(x_grid, function(x){
      bw <- x_grid[2] - x_grid[1]
      r <- range(x_grid)
      x_t1 <- sapply(z_m, function(z) f(x, z*h, p))
      a <- bin(z_g * z_i * x_t1, binwidth=bw, range=c(r[1], r[2]-bw))$count
      a / sum(a) 
    })
    t(mat)
  })
  SDP_Mat
}

# cubature



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
#'  matrix with the indices of h_grid giving the optimal h at each value x
#'  as the columns, with a column for each time.  
#'  V is a matrix of x_grid by x_grid, which is used to store the value 
#'  function at each point along the grid at each point in time.  
#'  The returned V gives the value matrix at the first (last) time.  
#' @export
find_dp_optim <- function(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, 
                          delta, reward=10){

 
  ## Initialize space for the matrices
  gridsize <- length(x_grid)
  HL <- length(h_grid)
  D <- matrix(NA, nrow=gridsize, ncol=OptTime)
  V <- rep(0,gridsize) # initialize BC,

  # give a fixed reward for having value larger than xT at the end. 
  V[ x_grid >= xT ] <- reward # a "scrap value" for x(T) >= xT

  # loop through time  
  for(time in 1:OptTime){
    # try all potential havest rates
    V1 <- sapply(1:HL, function(i){
      # Transition matrix times V gives dist in next time
      SDP_Mat[[i]] %*% V + 
      # then (add) harvested amount times discount
       profit(x_grid, h_grid[i]) * (1 - delta) 
## CHECKME I think discount is exp(-delta * (OptTime - time)) only in cts time
    })

    # find havest, h that gives the maximum value
    out <- sapply(1:gridsize, function(j){
      value <- max(V1[j,], na.rm = T) # each col is a diff h, max over these
      index <- which.max(V1[j,])  # store index so we can recover h's 
      c(value, index) # returns both profit value & index of optimal h.  
    })

    # Sets V[t+1] = max_h V[t] at each possible state value, x
    V <- out[1,]                        # The new value-to-go
    D[,OptTime-time+1] <- out[2,]       # The index positions
  }

  # Reed derives a const escapement policy saying to fish the pop down to
  # the largest population for which you shouldn't harvest: 
  ReedThreshold <- x_grid[max(which(D[,1] == 1))]

  # Format the output 
  list(D=D, V=V, S=ReedThreshold)
}



#' Forward simulate given the optimal havesting policy, D
#' @param f the growth function of the escapement population (x-h)
#'   should be a function of f(y, p), with parameters p
#' @param p the parameters of the growth function
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param Xo initial stock size
#' @param D the optimal solution indices on h_grid, 
#'  given for each possible state at each timestep
#' @param z_g a function which returns a random multiple for population growth 
#' @param z_m a function which returns the measurement uncertainty in the 
#'  assessment of stock size
#' @param z_i a function which returns a random number from a distribution 
#' for the implementation uncertainty in quotas 
#' @param interval is the years between updating the harvest quota
#' @returns a data frame with the time, fishstock, harvested amount,
#'  and what the stock would have been without that year's harvest. 
#' @export
ForwardSimulate <- function(f, pars, x_grid, h_grid, x0, D, z_g,
                            z_m, z_i, interval = 1){
  # initialize variables with initial conditions
  OptTime <- dim(D)[2]    # Stopping time
  x_h <- numeric(OptTime) # population dynamics with harvest
  x   <- numeric(OptTime) # What would happen with no havest
  h   <- numeric(OptTime) # optimal havest level
  x_h[1] <- x0  # initial values
  x[1]   <- x0  # intial values

  ## Simulate through time ##
  for(t in 1:(OptTime-1)){
    # Assess stock, with potential measurement error
    m_t <- x_h[t] * z_m()
    # Current state (is closest to which grid posititon) 
    St <- which.min(abs(x_grid - m_t)) 
    # Set harvest quota on update years
    if(t %% interval == 0)
      q_t <- h_grid[D[St, t + 1]] 
    else 
      q_t <- h[t]
    # Implement harvest/(effort) based on quota with noise 
    h[t + 1] <- q_t * z_i()
    # Noise in growth 
    z <- z_g() 
    # population grows
    x_h[t+1] <- z * f(x_h[t], h[t + 1], pars) # with havest
    x[t+1]   <- z * f(x[t], 0, pars) # havest-free dynamics
  }
  # formats output 
  data.frame(time=1:OptTime, fishstock=x_h, harvest=h, unharvested=x) 
}


