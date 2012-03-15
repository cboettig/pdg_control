# file determine_SDP_matrix.R 
# author Carl Boettiger <cboettig@gmail.com>
# date 2011-11-07
# licence BSD
# SDP solution adapted from SDP.m by Michael Bode
# 
# Contains three different functions for determining the 
# Stochastic Transition Matrix used in the Dynamic Programming solution
#


#########################################################################
# A function to generate the transition matrix used for the SDP routine #
#########################################################################

#' Determine the Stochastic Dynamic Programming matrix.
#' @param models: potential growth fns, must be function(x,h) with
#' pars hard-coded in.  
#' @param x_grid the discrete values allowed for the population size, x
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param sigma_g the variance of the population growth process
#' @return the transition matrix at each value of h in the grid. 
#' @details this analytical approach doesn't reliably support other 
#'  sources of variation.  The quality of the analytic approximations 
#'  (lognormal) can be tested. 
#' @export
determine_SDP_matrix <- function(models, x_grid, h_grid, p_grid, sigma_g){
  lapply(h_grid, function(h){
      transition <- 
      function(y,q,x,p){
        P <- lapply(models, function(f){
          fx <- f(x,h)
          x_to_y <- dlnorm(y, log(fx)-sigma_g^2/2, sigma_g) # logmean of 0 is log(1)
          x_to_y <- (fx <= 0 && y == 0 ) * 1 + (fx > 0) * x_to_y + 0 # transition all states to 0 if f(x) <= 0
        })
        p_to_q <- p * P[[1]] / (p * P[[1]] + (1-p) * P[[2]])
        P[[1]] * as.numeric( q == p_grid[ (p_grid > p_to_q)][1] )    # transition with prob p if q is correct, otherwise prob = 0 
      }
      out <-
      lapply(x_grid, function(y)
        sapply(p_grid, function(q)
          sapply(x_grid, function(x)
            sapply(p_grid, function(p)
              transition(y,q,x,p) ))))
  })
}




      nx <- length(x_grid)
      ny <- length(p_grid)
      Tmatrix <- matrix(NA, ncol = nx * ny, nrow = nx * ny)
      for(i in 1:nx)
        for(j in 1:ny)
          for(k in 1:nx)
            for(l in 1:ny)
              Tmatrix[i+(j-1)*ny, k+(l-1)*ny] <- 
                transition(x_grid[i], p_grid[j], x_grid[k], p_grid[l])
      Tmatrix


########################################################################
# A function to identify the dynamic optimum using backward iteration  #
########################################################################

#' Identify the dynamic optimum using backward iteration (dynamic programming)
#' @param Tmatrices a list of the stochastic transition matrix at each h value,
#' for each model i
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
find_dp_optim <- function(Tmatrices, x_grid, h_grid, OptTime, xT, profit, 
                          delta, reward=0){
  ## Initialize space for the matrices
  gridsize <- length(x_grid)
  HL <- length(h_grid)
  ML <- length(Tmatrices)
  D <- matrix(NA, nrow=gridsize, ncol=OptTime)
  V <- rep(0,gridsize) # initialize BC,

  # give a fixed reward for having value larger than xT at the end. 
  V[x_grid >= xT] <- reward # a "scrap value" for x(T) >= xT

  # loop through time  
  for(time in 1:OptTime){ 
    # try all potential havest rates
    V1 <- sapply(1:HL, function(h){
      sapply(1:ML, function(i)
      sigma[i, time] * Tmatrices[[i]][[h]] %*% V + 
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



