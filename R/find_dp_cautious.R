
#' Identify the best policy which harvests at P percent of optimum harvest
#' @param SDP_Mat the stochastic transition matrix at each h value
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param OptTime the stopping time 
#' @param xT the boundary condition population size at OptTime
#' @param c the cost/profit function, a function of harvested level
#' @param delta the exponential discounting rate
#' @param reward  the profit for finishing with >= Xt fish at the end 
#' (i.e. enforces the boundary condition)
#' @param P fraction of optimal harvest to use as the policy.  
#' @return list containing the matrices D and V.  D is an x_grid by OptTime
#'  matrix with the indices of h_grid giving the optimal h at each value x
#'  as the columns, with a column for each time.  
#'  V is a matrix of x_grid by x_grid, which is used to store the value 
#'  function at each point along the grid at each point in time.  
#'  The returned V gives the value matrix at the first (last) time. 
#' @export
find_dp_cautious <- function(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, 
                          delta, reward=0, P=1){

 
  ## Initialize space for the matrices
  gridsize <- length(x_grid)
  HL <- length(h_grid)
  D <- matrix(NA, nrow=gridsize, ncol=OptTime)
  V <- rep(0,gridsize) # initialize BC,

  ## give a fixed reward for having value larger than xT at the end. 
  V[x_grid >= xT] <- reward # a "scrap value" for x(T) >= xT

  ## loop through time  
  for(time in 1:OptTime){ 
    ## try all potential havest rates
    V1 <- sapply(1:HL, function(i){
      ## Transition matrix times V gives dist in next time
      SDP_Mat[[i]] %*% V + 
      ## then (add) harvested amount times discount
       profit(x_grid, h_grid[i]) * (1 - delta) 
    })

    ## find havest, h that gives the maximum value
    out <- sapply(1:gridsize, function(j){
      
      ## Optimal harvest solution is:
      best_index <- which.max(V1[j,])  

      ## Set policy to harvest level closest to fraction P of optimal harvest
      index <- which.min(abs(h_grid -  P * h_grid[best_index]))
      ## each col is a diff h, use the one corresponding to the policy chosen 
      value <- V1[j,index] 
      c(value, index) # returns both profit value & index of optimal h.  
    })

    ## Sets V[t+1] = max_h V[t] at each possible state value, x
    V <- out[1,]                        # The new value-to-go
    D[,OptTime-time+1] <- out[2,]       # The index positions
  }

  ## Reed derives a const escapement policy saying to fish the pop down to
  ## the largest population for which you shouldn't harvest: 
  ReedThreshold <- x_grid[max(which(D[,1] == 1))]

  ## Format the output 
  list(D=D, V=V, S=ReedThreshold)
}



