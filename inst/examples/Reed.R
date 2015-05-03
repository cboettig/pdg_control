# file Reed.R
# author Carl Boettiger, <cboettig@gmail.com>
# date 2011-11-02
# license BSD
# modified from Reed_SDP.m, by Michael Bode.  
# 
# Implements a numerical version of the SDP described in:
#   Reed, W.J., 1979. Optimal Escapement Levels in Stochastic
#   and Deterministic Harvesting Models. Journal of Environmental 
#   Economics and Management. 6: 350-363.
#
#           Fish population dynamics:
#           X_{t+1} = Z_n f(X_n) 


require(pdgControl)
require(reshape2)
require(ggplot2)

## consider defaults for these
# Define all parameters 
delta <- 0.05      # economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 10   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth

# load noise distributions  
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1


### Chose the state equation / population dynamics function
f <- BevHolt                # Select the state equation
pars <- c(1.5, 0.05)        # parameters for the state equation
xT <- 0                     # boundary conditions
x0 <- 10 
profit <- profit_harvest() 

# Set up the discrete grids
x_grid <- seq(0, 12, length = gridsize)  # population size
h_grid <- x_grid  # vector of havest levels, use same resolution as for stock

SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )

opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta)

ForwardSimulate(f,pars,x_grid,h_grid,x0, opt$D, z_g)

