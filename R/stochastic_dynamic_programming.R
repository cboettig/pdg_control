# file stochastic_dynamic_programming.R
# author Carl Boettiger <cboettig@gmail.com>
# date 2011-11-07
# licence BSD
# adapted from SDP.m by Michael Bode
# 
# Contains functions for stochastic dynamic programming 
# see Reed.R for example uses 

#' 


#' Define a profit function, price minus cost
#' @param p price of fish (Note, optimal will scrap xT if price is high enough!) 
#' @param c fishing extraction costs (per unit effort)
#' @return a function computing profit at effort intensity h_i over
#' possible stock values x_grid, profit(x_grid, h_i)
#' @export
profit_effort <- function(p = 1, c = 0.001){
#' @param x_grid is a the grid of state values (profit will evaluate at each of them)
#' @param h_i is the current harvesting *effort* (effort*stocksize = catch) 
#' @return the profits of fishing at intensity h_i given the stock value equals x_i
#' for each x_i in the grid.   
#' @details Due to the symmetry, you can actually compute the profit over a range 
#'  of harvest values, rather than a range of stock values, by simply swapping 
#'  x and h, i.e. give a vector of h values as x_grid, and a single stock size as h_i. 
  function(x_grid, h_i){ 
    harvest <- x_grid * h_i 
    sapply(harvest, function(x) max(0, p * x - c / x))
  }
}


#' Define a profit function, price minus cost
#' @param p price of fish (Note, optimal will scrap xT if price is high enough!) 
#' @param c fishing extraction costs (per unit effort)
#' @return a function computing profit at harvest intensity h_i over
#' possible stock values x_grid, profit(x_grid, h_i)
#' @export
profit_harvest  <- function(p = 1, c = 0.000){
#' @param x_grid is a the grid of state values (profit will evaluate at each of them)
#' @param h_i is total harvest level
#' @return the profits of harvesting at intensity h_i for each possible stock
#' size in x_grid.  
#' @details Due to the symmetry, you can actually compute the profit over a range 
#'  of harvest values, rather than a range of stock values, by simply swapping 
#'  x and h, i.e. give a vector of h values as x_grid, and a single stock size as h_i. 
  function(x_grid, h_i){
    sapply(x_grid, function(x_i){
      max(0, p * min(h_i, x_i) - c / x_i)
    })
  }
}




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
#' @return the transition matrix at each value of h in the grid. 
#' @details this analytical approach doesn't reliably support other 
#'  sources of variation.  The quality of the analytic approximations 
#'  (lognormal) can be tested. 
#' @export
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
        Prob <- dlnorm(ProportionalChance, 0, sigma_g) #logmean of 0 is log(1)
        # Store normalized probabilities in row
        SDP_matrix[i,] <- Prob/sum(Prob)
## dnorm(x_i, mu, sigma)  x_i+1,  
      }
    }
    SDP_matrix
  })
  SDP_Mat
}


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
#' @return the transition matrix at each value of h in the grid.  
#' @import snowfall
#' @import ggplot2
#' @import Hmisc
#' @export
SDP_by_simulation <- function(f, p, x_grid, h_grid, z_g, z_m, z_i, reps = 999){
  require(snowfall) # support parallelization of this
  sfExportAll()
  sfLibrary(ggplot2) # for the bin function 

  bw <- x_grid[2] - x_grid[1]
  lower <-x_grid[1]
  upper <- x_grid[length(x_grid)-1]

  SDP_Mat <- sfLapply(h_grid, function(h){ 
    mat <- sapply(x_grid, function(x){
      x_t1 <- replicate(reps, z_g() * z_i() * f(x, z_m() * h, p)) 
      a <- bin( x_t1, binwidth=bw, range=c(lower, upper))$count
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
#' @param reward  the profit for finishing with >= Xt fish at the end 
#' (i.e. enforces the boundary condition)
#' @return list containing the matrices D and V.  D is an x_grid by OptTime
#'  matrix with the indices of h_grid giving the optimal h at each value x
#'  as the columns, with a column for each time.  
#'  V is a matrix of x_grid by x_grid, which is used to store the value 
#'  function at each point along the grid at each point in time.  
#'  The returned V gives the value matrix at the first (last) time. 
#' @export
find_dp_optim <- function(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, 
                          delta, reward=0){

 
  ## Initialize space for the matrices
  gridsize <- length(x_grid)
  HL <- length(h_grid)
  D <- matrix(NA, nrow=gridsize, ncol=OptTime)
  V <- rep(0,gridsize) # initialize BC,

  # give a fixed reward for having value larger than xT at the end. 
  V[x_grid >= xT] <- reward # a "scrap value" for x(T) >= xT

  # loop through time  
  for(time in 1:OptTime){ 
    # try all potential havest rates
    V1 <- sapply(1:HL, function(i){
      # Transition matrix times V gives dist in next time
      SDP_Mat[[i]] %*% V + 
      # then (add) harvested amount times discount
       profit(x_grid, h_grid[i]) * (1 - delta) 
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
#' @return a data frame with the time, fishstock, harvested amount,
#'  and what the escapement ("unharvested"). 
#' @export
ForwardSimulate <- function(f, pars, x_grid, h_grid, x0, D, z_g,
                            z_m, z_i){
  # initialize variables with initial conditions
  OptTime <- dim(D)[2]    # Stopping time
  x_h <- numeric(OptTime) # population dynamics with harvest
  h <- numeric(OptTime) # optimal havest level
  x_h[1] <- x0  # initial values
  
  s <- x_h # also track escapement
  x <- x_h # What would happen with no havest
  
    
  ## Simulate through time ##
  for(t in 1:(OptTime-1)){
    # Assess stock, with potential measurement error
    m_t <- x_h[t] * z_m()
    # Current state (is closest to which grid posititon) 
    St <- which.min(abs(x_grid - m_t)) 
    q_t <- h_grid[D[St, (t + 1) ]] 
    # Implement harvest/(effort) based on quota with noise 
    h[t] <- q_t * z_i()
    # Noise in growth 
    z <- z_g() 
    # population grows
    x_h[t+1] <- z * f(x_h[t], h[t], pars) # with havest
    s[t+1]   <- x_h[t] - q_t # anticipated escapement
    x[t+1]   <- z * f(x[t], 0, pars) # havest-free dynamics
  }
  # formats output 
  data.frame(time=1:OptTime, fishstock=x_h, harvest=h, unharvested=x, escapement=s) 
}


########################################################################
# A function to identify the dynamic optimum using backward iteration  #
########################################################################
#' Optimum policy when changing harvest level has an adjustment cost
#' 
#' Identify the dynamic optimum using backward iteration (dynamic programming)
#' @param SDP_Mat the stochastic transition matrix at each h value
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param OptTime the stopping time 
#' @param xT the boundary condition population size at OptTime
#' @param c the cost/profit function, a function of harvested level
#' @param delta the exponential discounting rate
#' @param reward the profit for finishing with >= Xt fish at the end 
#' (i.e. enforces the boundary condition)
#' @param P the cost of adjusting a policy, proportional to the amount of change
#' @param penalty the kind of penalty applied: currently L1, L2, or assymetric
#' @return list containing the matrices D and V.  D is an x_grid by OptTime
#'  matrix with the indices of h_grid giving the optimal h at each value x
#'  as the columns, with a column for each time.  
#'  V is a matrix of x_grid by x_grid, which is used to store the value 
#'  function at each point along the grid at each point in time.  
#'  The returned V gives the value matrix at the first (last) time. 
#' @import expm
#' @export
optim_policy <- function(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, 
                          delta, reward=0, P=0, 
                          penalty=c("L1", "L2", "asymmetric", "fixed", "none")){
  penalty <- match.arg(penalty)
  ## Initialize space for the matrices
  gridsize <- length(x_grid)
  HL <- length(h_grid)
  D <- lapply(1:HL, function(i) matrix(NA, nrow=gridsize, ncol=OptTime))
  V <- sapply(1:HL, function(i){ 
    Vi <- rep(0,gridsize)  
    # give a fixed reward for having value larger than xT at the end. 
    Vi[x_grid >= xT] <- reward # a "scrap value" for x(T) >= xT
    Vi
  }) 

 
  # loop through time  
  for(time in 1:OptTime ){ 
  # loop over all possible values for last year's harvest level
   for(h_prev in 1:HL){
      # try all potential havest rates
      V1 <- sapply(1:HL, function(i){
        # cost of changing the policy from the previous year
        if(penalty=="L2")
          change_cost <- P * (h_grid[i] - h_grid[ h_prev])^2 
        else if(penalty=="L1")
          change_cost <- P * abs(h_grid[i] - h_grid[ h_prev]) 
        else if(penalty=="asymmetric")
          change_cost <- P * max(h_grid[h_prev]-h_grid[i], 0) 
        else if(penalty=="fixed")
          change_cost <- P
        else 
          0
        # Transition matrix times V gives dist in next time
        SDP_Mat[[i]]  %*% V[, h_prev] + 
        # then (add) harvested amount - policy cost times discount
         (profit(x_grid, h_grid[i]) - change_cost ) * (1 - delta)
      })
      # find havest, h that gives the maximum value
      out <- sapply(1:gridsize, function(j){
        value <- max(V1[j,], na.rm = T) # each col is a diff h, max over these
        index <- which.max(V1[j,])  # store index so we can recover h's 
        c(value, index) # returns both profit value & index of optimal h.  
      })
      # Sets V[t+1] = max_h V[t] at each possible state value, x
      V[, h_prev] <- out[1,]                        # The new value-to-go
      D[[h_prev]][,OptTime-time+1] <- out[2,]       # The index positions
    }
  }
  # Reed derives a const escapement policy saying to fish the pop down to
  # the largest population for which you shouldn't harvest: 
  # ReedThreshold <- x_grid[max(which(D[,1] == 1))]
  # Format the output 
  list(D=D, V=V)
}







#' Forward simulate given the optimal havesting policy when cost depends on prev harvest
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
#' @param alt_D an option to compute harvest with an alternate 
#' "optimal" solution as well, for comparison (indep of h_prev) 
#' @return a data frame with the time, fishstock, harvested amount,
#'  and what the escapement ("unharvested") is. Also the alternate
#'  fishstock and harvest dynamics, if calculated
#' @export
simulate_optim <- function(f, pars, x_grid, h_grid, x0, D, z_g,
                            z_m, z_i, alt_D = NULL, OptTime = dim(D[[1]])[2] ){
  # initialize variables with initial conditions
  x_h <- numeric(OptTime) # population dynamics with harvest
  h <- numeric(OptTime) # optimal havest level
  x_h[1] <- x0  # initial values
  h_prev <- 1 # assume no harvesting happening at start (index of h_grid)
  s <- x_h # also track escapement
  x <- x_h # What would happen with no havest

  # initialize alternative fishstock/harvest under alt_D management
  x_alt <- x_h
  h_alt <- h
    
  ## Simulate through time ##
  for(t in 1:(OptTime-1)){
    # Draw the noise random variables
    zg <- z_g() 
    zm <- z_m()
    zi <- z_i()
    # Assess stock, with potential measurement error
    m_t <- x_h[t] * zm
    # Current state (is closest to which grid posititon) 
    St <- which.min(abs(x_grid - m_t)) 
    # Set harvest quota on update years
    q_t <- h_grid[D[[h_prev]][St, (t + 1) ]] 
    # Implement harvest/(effort) based on quota with noise
    h[t] <- q_t * zi
    # Store last years quota
    h_prev <- which.min(abs(h_grid - q_t)) 
    # population grows
    x_h[t+1] <- zg * f(x_h[t], h[t], pars) # with havest
    s[t+1]   <- x_h[t] - q_t # anticipated escapement
    x[t+1]   <- zg * f(x[t], 0, pars) # havest-free dynamics

    ## Comparison solution
    if(!is.null(alt_D)){
      m_alt <- x_alt[t] * zm
      # Current state (is closest to which grid posititon) 
      St_alt <- which.min(abs(x_grid - m_alt)) 
      # Set harvest quota 
      h_alt[t] <- h_grid[alt_D[St_alt, (t + 1) ]] * zi
      # population grows
      x_alt[t+1] <- zg * f(x_alt[t], h_alt[t], pars) 
    }
}
  # formats output 
  data.frame(time=1:OptTime, fishstock=x_h, harvest=h, 
             unharvested=x, escapement=s, alternate=x_alt,
             harvest_alt = h_alt) 
}


