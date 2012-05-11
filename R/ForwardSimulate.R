## Note: See policy_costs.R for another simulation algorithm 


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
#' @param profit function, if known
#' @return a data frame with the time, fishstock, harvested amount,
#'  and what the escapement ("unharvested"). 
#' @export
ForwardSimulate <- function(f, pars, x_grid, h_grid, x0, D, z_g,
                            z_m=function(x) 1, z_i = function(x) 1, 
                            profit=NULL){
  # initialize variables with initial conditions
  OptTime <- dim(D)[2]    # Stopping time
  x_h <- numeric(OptTime) # population dynamics with harvest
  h <- numeric(OptTime) # optimal havest level
  x_h[1] <- x0  # initial values
  
  s <- x_h # also track escapement
  x <- x_h # What would happen with no havest
  p <- numeric(OptTime)
 
    
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
    p[t] <- profit(x_h[t], h[t])
}
  # formats output 
  data.frame(time = 1:OptTime, fishstock = x_h, harvest = h,
             unharvested = x, escapement = s, profit = p) 
}




