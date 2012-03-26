update_belief = function(x_t0, p_t0, x_t1){
  y1 = p_t0 * f1(x_t1, x_t0)
  y2 = (1-p_t0) * f2(x_t1, x_t0)
  P1 = y1 / (y1 + y2)
  if(is.na(P1) || x_t0 == 0)
    P1 = p_t0
  else{
    i = 1
    np = length(p_grid)
    while(p_grid[i] < P1 & i < np)
      i = i+1
    P1 = p_grid[i]  
 }
 P1
}




#' Forward simulate given the optimal havesting policy, D
#' @param f the true growth function of the escapement population (x-h)
#'   should be a function of f(y, p), with parameters p
#' @param pars the parameters of the growth function
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param p_grid the discrete values of belief allowed
#' @param x0 initial stock size
#' @param p0 initial belief
#' @param D the optimal solution indices on h_grid, 
#'  given for each possible state at each timestep
#' @param z_g a function which returns a random multiple for population growth 
#' for the implementation uncertainty in quotas 
#' @param update_belief a function of x[t], p[t], x[t+1]
#' @return a data frame with the time, fishstock, harvested amount,
#'  and what the escapement ("unharvested"). 
#' @export
active_adaptive_simulate <- function(f, pars, x_grid, h_grid, p_grid, x0, 
                                     p0, D, z_g, update_belief){
  # initialize variables with initial conditions
  OptTime <- dim(D)[2]    # Stopping time
  x_h <- numeric(OptTime) # population dynamics with harvest
  h <- numeric(OptTime)   # optimal havest level
  x_h[1] <- x0            # initial values
  p <- numeric(OptTime)   # belief
  p[1] <- p0              # initial belief
  s <- x_h                # also track escapement
  x <- x_h                # What would happen with no havest
  nx <- length(x_grid)
  np <- length(p_grid)
  getpolicy <- function(p,x, time) D[x + (p-1)*nx, time]
  ## Simulate through time ##
  for(t in 1:(OptTime-1)){
    # Current state (is closest to which grid posititon) 
    St <- which.min(abs(x_grid - x_h[t])) 
    h[t] <- h_grid[getpolicy(p[t], St, t)] 
    # Implement harvest/(effort) based on quota with noise 
    # Noise in growth 
    z <- z_g() 
    # population grows
    x_h[t+1] <- z * f(x_h[t], h[t], pars) # with havest
    s[t]     <- x_h[t] - h[t] # anticipated escapement
    x[t+1]   <- z * f(x[t], 0, pars) # havest-free dynamics
    p[t+1]   <- update_belief(x[t], p[t], x[t+1])
  }
  # formats output 
  data.frame(time = 1:OptTime, fishstock = x_h, harvest = h,
             unharvested = x, escapement = s) 
}


