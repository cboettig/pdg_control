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

#' Determine the Stochastic Dynamic Programming matrix when state space is 
#' multidimensional
#' @param f the growth function of the escapement population (x-h)
#'   should be a function of f(t, y, p), with parameters p, states y
#' @param p the parameters of the growth function
#' @param x_grid a data frame of discrete values allowed for each of the 
#'  state variables.  
#' @param h_grid the discrete values of harvest levels to optimize over
#' @param sigma_g the variance of the population growth process
#' @return the transition matrix at each value of h in the grid. 
#' @details this analytical approach doesn't reliably support other 
#'  sources of variation.  The quality of the analytic approximations 
#'  (lognormal) can be tested. 
#' @export
multid_SDP_matrix <- function(f, p, x_grid, h_grid, sigma_g){


  gridsizes <- c(a=5, s=6, x=7 ,h=4)

  sigma_g <- 0.2
  a_grid <- seq(1, 2, length=gridsizes["a"])
  s_grid <- seq(.0, .1, length=gridsizes["s"])  
  x_grid <- seq(0, 1, length=gridsizes["x"]) 
  h_grid <- seq(0, .5, length=gridsizes["h"])


  f <- function(x,a,h,p){
    max(0, (x - h) * a / (1 + (x - h) * p[1]))
  }
  p = c(1)

## For each value of x_i, find the probability of going to x_j (given a_i and s_i)
# P( x_i -> x_j; | a_i, s_i)
    px <- function(a_i, s_i, x_i)
      sapply(x_grid, function(x_j){
        sum(
          sapply(a_grid, function(a){
            dlnorm(x_j, log(f(x_i, a, h, p)) - sigma_g ^ 2 / 2, sigma_g) *
            dnorm(a, a_i, s_i)
          })
        )
      })
  x_
      px(a_i, s_i, x_t_minus_1)
  
  # a_i -> a_j | x_i, s_i
P_to_aj_given <- 

  
  sapply(a_grid, function(a_j){ 
       dlnorm( log(f(x_tminus1_i, a_j, h,p) - sigma_g ^ 2 / 2, sigma_g)) * dnorm(a_j, a_i, s_i)
      })

P_to_sj_given <- 


dat <-
sapply(h_grid, function(h){
  sapply(a_grid, function(a_hat){
    sapply(s_grid, function(sigma_a){
      sapply(x_grid, function(x){
        sapply(a_grid, function(a){
          sapply(x_grid, function(y){
            dlnorm(y, log(f(x, a, h, p))-sigma_g^2/2, sigma_g) *
            dnorm(a, a_hat, sigma_a)
          })
        })
      })
    })
  })
})


get_column <- function(a,s,x,h, dat, gridsizes){
  stride_a <- gridsize["x"] * gridsizes["a"] * gridsizes["x"] * gridsizes["s"]
  stride_s <- gridsize["x"] * gridsizes["a"] * gridsizes["x"]
  stride_x <- gridsize["x"] * gridsizes["a"]
  startpt <- 
  1 + (a-1) * stride_a +
  1 + (s - 1) * stride_s +
  1 + (x - 1) * stride_x
  endpt <- startpt + gridsizes["a"]
  dat[startpt:endpt, h]
}


get_SDP_matrix <- function(){
  n <- dim(dat)[1]  
  SDP_Mat <- matrix(NA, ncol=n, nrow=n)
        a1s1x1 
  a1s1x1 get_column(1,1,1, 1, dat, gridsizes)
}


# Format of h_a_s_x_A:
#    h1,         h2      h3 
#  a1,x1,s1,a1
#  a2,x1,s1
#  a3,x1,s1
#  a1,x2,s1
#  a2,x2,s1
#  a3 x2,s1
#  a1 x3,s1
#   ...
#  a1,x1,s2

find_dp_optim <- function(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, 
                          delta, reward=0){
  ## Initialize space for the matrices
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



