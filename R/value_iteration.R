
#' Identify the dynamic optimum using backward iteration (dynamic programming)
#' @param SDP_Mat the stochastic transition matrix at each h value
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param OptTime the stopping time 
#' @param xT the boundary condition population size at OptTime
#' @param c the cost/profit function, a function of harvested level
#' @param delta the discounting rate (1-delta)
#' @param epsilon value iteration tolerance

#' (i.e. enforces the boundary condition)
#' @return list containing the matrices D and V.  D is an x_grid by OptTime
#'  matrix with the indices of h_grid giving the optimal h at each value x
#'  as the columns, with a column for each time.  
#'  V is a matrix of x_grid by x_grid, which is used to store the value 
#'  function at each point along the grid at each point in time.  
#'  The returned V gives the value matrix at the first (last) time. 
#' @export
value_iteration <- function(SDP_Mat, x_grid, h_grid, OptTime=100, xT, profit, 
                            delta, epsilon = 1e-4){
  
  
  ## Initialize space for the matrices
  gridsize <- length(x_grid)
  HL <- length(h_grid)
  D <- matrix(NA, nrow=gridsize, ncol=1)
  V <- rep(0,gridsize) # initialize BC,
  Vp <- rep(10, gridsize) # initialize previous value (not equal to current value)
  
  # loop through time  
  #for(time in 1:OptTime){
  time <- 1
  while(max(abs(Vp - V)) >= (delta * epsilon / (2 - 2 * delta)) && time < OptTime){
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
    Vp <- V
    V <- out[1,]       # The new value-to-go
    D <- out[2,]       # The index positions
    
    
    time <- time+1
  }
  # check for convergence in V
  
  
  # Reed derives a const escapement policy saying to fish the pop down to
  # the largest population for which you shouldn't harvest: 
  ReedThreshold <- x_grid[max(which(D == 1))]
  
  # Format the output 
  list(D=D, V=V, S=ReedThreshold)
}


