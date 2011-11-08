# file Reed.R
# author Carl Boettiger, <cboettig@gmail.com>
# date 2011-11-02
# modified from Reed_SDP.m, by Michael Bode.  
# 
# Implements a numerical version of the SDP described in:
#   Reed, W.J., 1979. Optimal Escapement Levels in Stochastic
#   and Deterministic Harvesting Models. Journal of Environmental 
#   Economics and Management. 6: 350-363.
# 
# 
# Fish population dynamics:
# X_{t+1} = Z_n f(X_n) 
# 

#pars <- c(2,4)    # Beverton-Holt/f(x) pars, A, B
#K <- (pars[1]-1)/pars[2]   # Unharvested deterministic equib pop


rm(list=ls())   # Start wtih clean workspace 
require(ggplot2) # nicer plotting package

set.seed(1)
# load my script defining the SDP functions:
# determine_SDP_matrix, find_dp_optim, & ForwardSimulate
source("stochastic_dynamic_programming.R")

# Define all parameters 
delta <- 0.1      # economic discounting rate
OptTime <- 25     # stopping time
sigma <- 0.002      # Noise process
gridsize <- 100   # gridsize (discretized population)

# Chose the state equation / population dynamics function
f <- RickerAllee
pars <- c(1, 100, 30)
K <- 100

# define a profit function, price minus cost
profit <- function(h){
  p <- 1
  c <- 0.001 # higher extraction costs result in less fishing 
  # sapply: support for vector-valued h
  sapply(h, function(h) max(0,p*h - c/h))
}

# Set up the grid 
x_grid <- seq(0, 2*K-2, length=gridsize)  # population size
h_grid <- x_grid  # vector of havest levels, use same res as stock

# Calculate the transition matrix 
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma)

# Find the optimum by dynamic programming 
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, 30, profit, delta)

## Plot the results of a single run, against unharvested version  
out <- ForwardSimulate(f, pars, x_grid, h_grid, sigma, K/2, opt$D)
dat <- melt(out, id="time")
ggplot(dat, aes(time, value, color=variable)) + geom_line()


#######################################################################
# Now we'll simulate this process many times under this optimal havest#
#######################################################################
sims <- lapply(1:100, function(i){
# simulate the optimal routine on a stoch realization of growth dynamics
    sim <- ForwardSimulate(f, pars, x_grid, h_grid, sigma, K/2, opt$D)
    list(fishstock=sim$fishstock, unharvested=sim$unharvested)
})

# some reformatting
fished <- sapply(sims, function(x) x$fishstock)
optimal_havest <- melt(data.frame(year = 1:OptTime,fished), id="year")

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
unharvested <- melt(data.frame(year=1:OptTime, unfished), id="year")

# without fishing, how many populations have crashed
crashed <- sum(unfished[OptTime-1,]<pars[3])

p2 <- ggplot(unharvested,aes(year, value)) + 
  geom_line(aes(group=variable), col="gray") + 
  geom_line(aes(year, rowMeans(unfished)))  # Mean path

p2 <- p2 + opts(title=sprintf("Unfished dynamics, %d populations crash",
  crashed))
ggsave("unfished.png")







#ReedThreshold <- n[ sum(D[,1]==1) ]

# check out the deterministic dynamics
#x <- numeric(OptTime)
#x[1] <- 50
#for(t in 1:(OptTime-1))
#  x[t+1] <- f(x[t],pars)

