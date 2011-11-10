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
rm(list=ls())   # Start wtih clean workspace 
require(ggplot2) # nicer plotting package

set.seed(1)
# load my script defining the SDP functions:
# determine_SDP_matrix, find_dp_optim, & ForwardSimulate
source("stochastic_dynamic_programming.R")

# Define all parameters 
delta <- 0.1      # economic discounting rate
OptTime <- 25     # stopping time
sigma <- 0.2      # Noise process
gridsize <- 100   # gridsize (discretized population)

# Chose the state equation / population dynamics function
# f <- BevHolt #pars <- c(2,4) #K <- (pars[1]-1)/pars[2]
# f <- RickerAllee # pars <- c(1, 100, 30) # K <- 100
f <- Myer
pars <- c(1,2,6) 
p <- pars # shorthand 
K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
# Boundary value conditions
x0 <- K
xT <- p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 # allee threshold
e_star <- (p[1]*sqrt(p[3])-2)/2 ## Bifurcation point 


# define a profit function, price minus cost
profit <- function(h){
  p <- 1     # higher profits mean more fishing. (scraps xT if too high!) 
  c <- 0.001 # higher extraction costs result in less fishing 
  # sapply: support for vector-valued h
  sapply(h, function(h) max(0,p*h - c/h))
}

# Set up the grid 
x_grid <- seq(0, 2*K, length=gridsize)  # population size
#h_grid <- x_grid  # vector of havest levels, use same res as stock
h_grid <- seq(0, 2, length=gridsize) # Myers model based on harvesting effort 

# Calculate the transition matrix 
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma)

# Find the optimum by dynamic programming 
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta)


## Example plot the results of a single run, against unharvested version  
out <- ForwardSimulate(f, pars, x_grid, h_grid, sigma, x0, opt$D)
dat <- melt(out, id="time")
p0 <- ggplot(dat, aes(time, value, color=variable)) + geom_line() +  geom_abline(intercept=opt$S, slope=0, col="black") #+
p0 <- p0 + geom_abline(intercept=e_star, slope=0, col="green", lty = 2)

#  geom_line(data=subset(dat, variable=="harvest"), 
#  aes(time, value+opt$S), col="black")

#######################################################################
# Now we'll simulate this process many times under this optimal havest#
#######################################################################
sims <- lapply(1:100, function(i){
# simulate the optimal routine on a stoch realization of growth dynamics
    sim <- ForwardSimulate(f, pars, x_grid, h_grid, sigma, x0, opt$D)
    list(fishstock=sim$fishstock, unharvested=sim$unharvested)
})

# some reformatting
fished <- sapply(sims, function(x) x$fishstock)
optimal_havest <- melt(data.frame(year = 1:OptTime,fished), id="year")

# After optimal fishing, how many populations have crashed 
optimal_crashed <- sum(fished[OptTime-1,]<=xT)

# Assemble the plot
p1 <- ggplot(optimal_havest, aes(year, value)) + 
  geom_line(aes(group = variable), col = "gray") + 
  geom_line(aes(year, rowMeans(fished)))  + # Mean path
  geom_abline(intercept=opt$S, slope=0)
p1 <- p1 + opts(title=sprintf("Optimally Havested, %d populations crash",
  optimal_crashed))

# reformatting for the unhavested dynamics
unfished <- sapply(sims, function(x) x$unharvested)
unharvested <- melt(data.frame(year=1:OptTime, unfished), id="year")

# without fishing, how many populations have crashed
crashed <- sum(unfished[OptTime-1,]<=xT)
p2 <- ggplot(unharvested,aes(year, value)) + 
  geom_line(aes(group=variable), col="gray") + 
  geom_line(aes(year, rowMeans(unfished)))  # Mean path
p2 <- p2 + opts(title=sprintf("Unfished dynamics, %d populations crash",
  crashed))



#ggsave("samplerun.png", plot=p0)
#ggsave("fished.png", plot=p1)
#ggsave("unfished.png", plot=p2)








# check out the deterministic dynamics
#x <- numeric(OptTime)
#x[1] <- K
#for(t in 1:(OptTime-1))
#  x[t+1] <- f(x[t],0,pars)

