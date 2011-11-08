# file Reed.R
# author Carl Boettiger, <cboettig@gmail.com>
# date 2011-11-02
# modified from SDP.m, by Michael Bode.  
# 
# Implements a numerical version of the SDP described in:
#   Reed, W.J., 1979. Optimal Escapement Levels in Stochastic
#   and Deterministic Harvesting Models. Journal of Environmental 
#   Economics and Management. 6: 350-363.
# 
# 
# Fish population dynamics:
# X_{t+1} = Z_n f(X_n) 
# f(x) = 
# 
# 


require(ggplot2) # nicer plotting package
rm(list=ls())   # Start wtih clean workspace 

# load my script defining the SDP functions:
# determine_SDP_matrix, find_dp_optim, & ForwardSimulate
source("stochastic_dynamic_programming.R")


# temporal variation in f is possible, but increases the memory required, 
# though not the time for the optimization(?)

########## Define our population dynamics / state equation  ################ 
#' Harvested Beverton Holt growth model
#' @param x fish population that reproduces (usually x-h)
#' @param p parameters of the growth function, c(A, B), where
#'  A is the maximum growth rate and B the half-maximum in Beverton-Holt.  
#' @returns population next year
BevHolt <- function(x, p){
  A <- p[1] 
  B <- p[2] 
  sapply(x, function(x){ # use sapply so fn accepts vector-valued x
    x <- max(0, x)
    max(0, A*x/(1+B*x))
  })
}

#' Discrete-time model with an allee effect for alpha > 1
#' @param x the current population level
#' @param p vector of parameters c(r, alpha, K) 
#' @returns the population level in the next timestep
#' @details A Beverton-Holt style model with Allee effect
#' try with pars = c(1,2,100)
AlleeModel <- function(x, p){
   p[1] * x ^ p[2] / (1 + x ^ p[2] / p[3]) 
}

#' @param x the current population level
#' @param p a vector of parameters c(r, K, C) 
#' @returns the population level in the next timestep
RickerAllee <- function(x, p){
    x * exp(p[1] * (1 - x / p[2]) * (x - p[3]) / p[2] ) 
}



# Define all parameters 
delta <- 0.1      # economic discounting rate
OptTime <- 25     # stopping time
sigma <- 0.2      # Noise process
gridsize <- 100   # gridsize (discretized population)

# Chose the state equation / population dynamics function
f <- RickerAllee
pars <- c(1, 100, 30)
K <- 100

#pars <- c(2,4)    # Beverton-Holt/f(x) pars, A, B
#K <- (pars[1]-1)/pars[2]   # Unharvested deterministic equib pop

# define a profit function, price minus cost
profit <- function(h){
  p <- 1
  c <- 0.001 # higher extraction costs result in less fishing 
  # sapply: support for vector-valued h
  sapply(h, function(h) max(0,p*h - c/h))
}



# Set up the grid 
x_grid <- seq(0, 2*K, length=gridsize)  # population size
h_grid <- x_grid  # vector of havest levels, use same res as stock

# Calculate the transition matrix 
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma)
# Find the optimum by dynamic programming 
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, 0, profit, delta)



sims <- lapply(1:100, function(i){
# simulate the optimal routine on a stoch realization of growth dynamics
    sim <- ForwardSimulate(f, pars, x_grid, h_grid, sigma, K/2, opt$D)
    list(fishstock=sim$fishstock, unharvested=sim$unharvested)
})

# some reformatting
fished <- sapply(sims, function(x) x$fishstock)
dat <- data.frame(year = 1:OptTime,out)
optimal_havest <- melt(dat, id="year")

# After optimal fishing, how many populations have crashed 
optimal_crashed <- sum(fished[OptTime-1,]<pars[3])

# Assemble the plot
p1 <- ggplot(optimal_havest, aes(year, value)) + 
  geom_line(aes(group = variable), col = "gray") + 
  geom_line(aes(year, rowMeans(fished)))  # Mean path
# modify plot appearance
p1 <- p1 + opts(title=sprintf("Optimally Havested, %d populations crash",
  optimal_crashed))
ggsave("fished.png")

# reformatting for the unhavested dynamics
unfished <- sapply(sims, function(x) x$unharvested)
unharvested <- melt(data.frame(year=1:OptTime, free), id="year")

# without fishing, how many populations have crashed
crashed <- sum(unfished[OptTime-1,]<pars[3])

p2 <- ggplot(unharvested,aes(year, value)) + 
  geom_line(aes(group=variable), col="gray") + 
  geom_line(aes(year, rowMeans(unfished)))  # Mean path

p2 <- p2 + opts(title=sprintf("Unfished dynamics, %d populations crash",
  crashed))
ggsave("unfished.png")






## Plot the results of a single run, against unharvested version  
#out <- ForwardSimulate(f, pars, x_grid, h_grid, sigma, K/2, opt$D)
#ggplot(out) + geom_line(aes(time, fishstock)) + geom_line(aes(time, unharvested), col=I("green"))

#ReedThreshold <- n[ sum(D[,1]==1) ]

# check out the deterministic dynamics
#x <- numeric(OptTime)
#x[1] <- 50
#for(t in 1:(OptTime-1))
#  x[t+1] <- f(x[t],pars)

