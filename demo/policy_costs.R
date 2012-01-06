# file Reed.R
# author Carl Boettiger, <cboettig@gmail.com>
# date 2011-11-02
# license BSD
rm(list=ls())   # Start wtih clean workspace 
require(pdgControl)

# shouldn't be necessary - these are imported!
require(ggplot2)
require(Hmisc)
require(expm)

## consider defaults for these
# Define all parameters 
delta <- 0.04      # economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 40   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
interval <- 1     # period of updating the stock assessment

# Define noise distributions  
z_g <- function() rlnorm(1,  0, sigma_g) # growth noise
z_m <- function() 1 # no measurement (stock assessment) noise
z_i <- function() 1 # no implementation (quota) noise

## Chose the state equation / population dynamics function
f <- BevHolt              # Select the state equation
pars <- c(2,4)            # parameters for the state equation
K <- (pars[1]-1)/pars[2]  # Carrying capacity 
xT <- 0                   # boundary conditions
scrap_value <- 0          # reward profit offered for finishing >= xT stock
e_star <- 0               # model's bifurcation point (just for reference)
control <- "harvest"      # control variable is total harvest, h = e * x
x0 <- K/2                 # initial condition
# use a harvest-based profit function with default parameters
profit <- profit_harvest() # functions defined in stochastic_dynamic_programming.R

# Set up the discrete grids
x_grid <- seq(0, 2 * K, length = gridsize)  # population size
h_grid <- x_grid  # vector of havest levels, use same resolution as for stock

SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g)

## solution with policy cost P: 
opt <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, .25*K, 
                    profit, delta, reward=100, P=2*9.3, penalty="asym")

# Calculate the Reed optimum (e.g. no cost to policy adjustment):
reed <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, .25*K, 
                      profit, delta, reward=100, interval=interval)

sims <- lapply(1:100, function(i)
  simulate_optim(f, pars, x_grid, h_grid, x0, opt$D, z_g, z_m, z_i, reed$D)
)

## Prove it adjusts less
## add a function to compute total profit when charged for adjustments, prove Reed is suboptimal

## sanity check when P=0
#identical(sims[[1]][[2]], sims[[1]][[6]])



## Reshape and summarize data ###
dat <- melt(sims, id="time") # reshapes the data matrix to "long" form
## Show dynamics of a single replicate 
ex <- sample(1:100,1) # a random replicate
example <- subset(dat, variable %in% c("fishstock", "alternate","harvest", "harvest_alt") & L1 == ex)
example[[2]] <- as.factor(as.character(example[[2]]))
p0 <- ggplot(example) +
      geom_line(aes(time, value, color = variable), position="jitter") 
p0 <- p0 + geom_abline(intercept = reed$S, slope = 0, col = "darkred") 
print(p0)



p1 <- plot_replicates(sims)

