#' @export
setmodel <- function(f, pars){
  function(x_t1, x_t0){
    mu = f(x_t0, 0, pars)
    (mu <= 0) * (x_t1 == 0) +
    (mu > 0) * dlnorm(x_t1, log(mu), sigma_g)
  }
}

# internal function
bin <- function(v, v_grid){
  v_grid[which.min(abs(v_grid - v))]
#  i = 1
#  n = length(v_grid)
#  while(v_grid[i] < v & i < n)
#    i = i+1
#  v_grid[i]  
}

#' @export
update_belief = function(f1,f2){
  function(x_t0, p_t0, x_t1, p_grid){
    y1 = p_t0 * f1(x_t1, x_t0)
    y2 = (1-p_t0) * f2(x_t1, x_t0)
    P1 = y1 / (y1 + y2)
    if(is.na(P1) || x_t0 == 0)
      P1 = p_t0
    else{
      P1 = bin(P1, p_grid)
   }
   P1
  }
}


# internal function
transition = function(p_grid, f1, f2){
  function(x_t0, p_t0, x_t1, p_t1){
    np = length(p_grid)
    y1 = p_t0 * f1(x_t1, x_t0)
    y2 = (1-p_t0) * f2(x_t1, x_t0)
    P1 = y1 / (y1 + y2)
    if(is.na(P1) || x_t0 == 0)
      P1 = p_t0
    else
      P1 = bin(P1, p_grid)
   (y1+y2) * ( p_t1 == P1)
  }
}

#' @export
model_uncertainty <- function(f1, f2, x_grid, p_grid, h_grid){
  f <- transition(p_grid, f1, f2)
  lapply(h_grid, function(h){
    x_minus_h <- (x_grid-h) * as.integer( (x_grid-h)>0 )
    d = expand.grid(x_t0 = x_minus_h, p_t0 = p_grid, x_t1 = x_grid, p_t1 = p_grid)
    M = matrix(mapply(f, d$x_t0, d$p_t0, d$x_t1, d$p_t1), nrow = length(p_grid) * length(x_grid) )
    for(i in 1:dim(M)[1]) # normalize
      M[i,] = M[i,]/sum(M[i,])
    M
  })
}


#' @export
dp_optim <- function(M, x_grid, h_grid, OptTime, xT, profit, 
                          delta, reward=0, p_grid){
  nx <- length(x_grid)
  np <- length(p_grid)
  gridsize <- np * nx 
  HL <- length(h_grid)
  D <- matrix(NA, nrow=gridsize, ncol=OptTime)
  V <- rep(0,gridsize) # initialize BC,

  profit.grid <- function(x_grid, h_i)
    expand.grid(profit(x_grid, h_i), p_grid)[[1]]

  # give a fixed reward for having value larger than xT at the end. 
  indices_fixed_x <- function(x) (1:np-1)*nx + x
  V[sapply(x_grid[x_grid>xT], function(x) indices_fixed_x(x))] <- reward

  # loop through time  
  for(time in 1:OptTime){ 
    # try all potential havest rates
    V1 <- sapply(1:HL, function(i){
      # Transition matrix times V gives dist in next time
      M[[i]] %*% V + 
      # then (add) harvested amount times discount
       profit.grid(x_grid, h_grid[i]) * (1 - delta) 
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
  # Format the output 
  list(D=D, V=V)
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
  getpolicy <- function(p,x, time){ 
    i_x = which.min(abs(x_grid-x))
    i_p = which.min(abs(p_grid-p))
    D[i_x + (i_p-1)*nx, time]
  }
  ## Simulate through time ##
  for(t in 1:(OptTime-1)){
    # Current state (is closest to which grid posititon) 
    h[t] <- h_grid[getpolicy(p[t], x_h[t], t)] 
    # Implement harvest/(effort) based on quota with noise 
    # Noise in growth 
    z <- z_g() 
    # population grows
    x_h[t+1] <- z * f(x_h[t], h[t], pars) # with havest
    s[t]     <- x_h[t] - h[t] # anticipated escapement
    x[t+1]   <- z * f(x[t], 0, pars) # havest-free dynamics
    p[t+1]   <- update_belief(x[t], p[t], x[t+1], p_grid)
  }
  # formats output 
  data.frame(time = 1:OptTime, fishstock = x_h, harvest = h,
             unharvested = x, escapement = s, belief=p) 
}



