## Modified optimization and simulation techniques that allow for a cost to policy 

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

  profit_matrix <- sapply(h_grid, function(h) profit(x_grid, h))

  # loop through time  
  for(time in 1:OptTime ){ 
  # loop over all possible values for last year's harvest level
   for(h_prev in 1:HL){
      # try all potential havest rates
      V1 <- sapply(1:HL, function(i){
        # cost of changing the policy from the previous year
        sales <- profit_matrix[, i]
        change_cost <- 
        switch(penalty,
          L2 =  sales - P * (h_grid[i] - h_grid[ h_prev]) ^ 2,
          L1 = sales - P * abs(h_grid[i] - h_grid[ h_prev]),
          asymmetric = sales - P * max(h_grid[h_prev] - h_grid[i], 0),
          fixed = sales - P,
          none = 0)

        # Transition matrix times V gives dist in next time
        SDP_Mat[[i]]  %*% V[, h_prev] + change_cost * {1 - delta} 
        # then (add) harvested amount - policy cost times discount
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
#'
#' This simulates a process where there is a cost to changing the policy.  
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
