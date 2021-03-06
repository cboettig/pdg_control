# file Reed.R
# author Carl Boettiger, <cboettig@gmail.com>
# date 2011-11-02
# license BSD
# modified from Reed_SDP.m, by Michael Bode.  
# 
# Implements a numerical version of the SDP described in:
# 
#   Sethi, G., Costello, C., Fisher, A., Hanemann, M., & Karp, L. (2005). 
#   Fishery management under multiple uncertainty. Journal of Environmental
#   Economics and Management, 50(2), 300-318. doi:10.1016/j.jeem.2004.11.005
#
#   Reed, W.J., 1979. Optimal Escapement Levels in Stochastic
#   and Deterministic Harvesting Models. Journal of Environmental 
#   Economics and Management. 6: 350-363.
#
# 
#           Fish population dynamics:
#           X_{t+1} = Z_n f(X_n) 


rm(list=ls())   # Start wtih clean workspace 
require(pdgControl)
require(reshape2)
require(ggplot2)


## consider defaults for these
# Define all parameters 
delta <- 0.1      # economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota

# load noise distributions  
source("noise_dists.R")

## Chose the state equation / population dynamics function
f <- BevHolt                # Select the state equation
pars <- c(2, 4)             # parameters for the state equation
K <- (pars[1] - 1)/pars[2]  # Carrying capacity 
xT <- 0                     # boundary conditions
e_star <- 0                 # model's bifurcation point (just for reference)
control = "harvest"         # control variable is total harvest, h = e * x

## An alternative state equation, with allee effect: (uncomment to select)
#f <- Myer_harvest
#pars <- c(1, 2, 6) 
#p <- pars # shorthand 
#K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
#xT <- p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 # allee threshold
#e_star <- (p[1] * sqrt(p[3]) - 2) / 2 ## Bifurcation point 
#control <- "harvest"          # control variable is harvest effort, e = h / x (for price eqn)

# initial condition near equib size (note the stochastic deflation of mean)
x0 <- K - sigma_g ^ 2 / 2 

# use a harvest-based profit function with default parameters
##profit <- profit_harvest(c = 1) # we'll use two different versions and compare, see below

# Set up the discrete grids
x_grid <- seq(0, 2 * K, length = gridsize)  # population size
h_grid <- x_grid  # vector of havest levels, use same resolution as for stock
#h_grid <- seq(0, 2, length=gridsize) 


#######################################################################
# Calculate the transition matrix (with noise in growth only)         #
#######################################################################
## Fast & decent approximation (lognormal growth noise only)
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g)

## Most accurate way: use integral (lognormal growth noise only)
## didn't work to integrate all noise forms
#SDP_Mat <- integrate_SDP_matrix(f, pars, x_grid, h_grid, sigma_g)

## calculate the transition matrix by simulation, generic to types of noise
#require(snowfall) # use parallelization since this can be slow
#sfInit(parallel=TRUE, cpu=4)
#SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps=999)

#######################################################################
# Find the optimum by dynamic programming                             #
#######################################################################
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit_harvest(c=1), delta, reward=10)

## Compare to case with profit function having little extraction cost

cheap <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit_harvest(c=.0001), delta, reward=10)


opt_sim <-  ForwardSimulate(f, pars, x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
cheap_sim <-  ForwardSimulate(f, pars, x_grid, h_grid, x0, cheap$D, z_g, z_m, z_i)


# What if parameter estimation is inaccurate? e.g.:
#pars[1] <- pars[1] * 0.95


#######################################################################
# Now we'll simulate this process many times under this optimal havest#
#######################################################################
sims <- lapply(1:100, function(i){
# simulate the optimal routine on a stoch realization of growth dynamics
  ForwardSimulate(f, pars, x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})

sims2 <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0, cheap$D, z_g, z_m, z_i)
})
#######################################################################
# Plot the results                                                    #
#######################################################################

## Reshape and summarize data ###
dat <- melt(sims2, id="time") # reshapes the data matrix to "long" form
## Show dynamics of a single replicate 
ex <- sample(1:100,1) # a random replicate
example <- subset(dat, variable %in% c("fishstock", "alternate","harvest", "harvest_alt") & L1 == ex)
example[[2]] <- as.factor(as.character(example[[2]]))
p0 <- ggplot(example) +
      geom_line(aes(time, value, color = variable), position="jitter") 
p0 <- p0 + geom_abline(intercept = cheap$S, slope = 0, col = "darkred") 
print(p0)



p1 <- plot_replicates(sims)

