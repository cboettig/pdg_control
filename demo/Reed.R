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
reward <- 0       # bonus for satisfying the boundary condition

# load noise distributions  
source("noise_dists.R")

### Chose the state equation / population dynamics function
#f <- BevHolt                # Select the state equation
#pars <- c(2, 4)             # parameters for the state equation
#K <- (pars[1] - 1)/pars[2]  # Carrying capacity 
#xT <- 0                     # boundary conditions
#e_star <- 0                 # model's bifurcation point (just for reference)
#control = "harvest"         # control variable is total harvest, h = e * x
#price <- 1
#cost <- K^2/2^2
#
## An alternative state equation, with allee effect: (uncomment to select)
f <- Myer_harvest
pars <- c(1, 2, 6) 
p <- pars # shorthand 
K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
xT <- p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 # allee threshold
e_star <- (p[1] * sqrt(p[3]) - 2) / 2 ## Bifurcation point 
control <- "harvest"          # control variable is harvest effort, e = h / x (for price eqn)
price <- 1
cost <- 0.00001

# initial condition near equib size (note the stochastic deflation of mean)
x0 <- K - sigma_g ^ 2 / 2 

# use a harvest-based profit function with default parameters
profit <- profit_harvest(p=price, c = cost) 

# Set up the discrete grids
x_grid <- seq(0, 2 * K, length = gridsize)  # population size
h_grid <- x_grid  # vector of havest levels, use same resolution as for stock
#h_grid <- seq(0.01, K, length=gridsize) 


#######################################################################
# Calculate the transition matrix (with noise in growth only)         #
#######################################################################
## Fast & decent approximation (lognormal growth noise only) # over or understimating noise does what you think
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )

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
                     profit, delta, reward)


# What if parameter estimation is inaccurate? e.g.:
#pars[1] <- pars[1] * 0.95


#######################################################################
# Now we'll simulate this process many times under this optimal havest#
#######################################################################
sims <- lapply(1:100, function(i){
# simulate the optimal routine on a stoch realization of growth dynamics
  ForwardSimulate(f, pars, x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})

#######################################################################
# Plot the results                                                    #
#######################################################################

# Make data tidy
dat <- melt(sims, id=names(sims[[1]]))  
# Faster if we use data.table objects
require(data.table)
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice

## Compute some further descriptors
crashed <- dt[time==OptTime, fishstock == 0, by=reps]
rewarded <- dt[time==OptTime, fishstock > xT, by=reps]

dt <- data.table(dt, id=1:5000)
profits <- dt[, profit(fishstock, harvest), by=id]

### And add this information to the data.table
setkey(dt, id)
setkey(profits, id)
dt <- dt[profits]
setnames(dt, "V1", "profits")
setkey(dt, reps)

## Compute total profit
total_profit <- dt[,sum(profits), by=reps]
## Add the boundary reward into the profit total?
total_profit <- total_profit + rewarded$V1 * reward 

## Add these three columns to the data.table
setkey(total_profit, reps)
setkey(crashed, reps)
setkey(rewarded, reps)
dt <- dt[total_profit]
dt <- dt[crashed]
dt <- dt[rewarded]
setnames(dt, c("V1", "V1.1", "V1.2"), c("total.profit", "crashed", "rewarded"))


p1 <- ggplot(dt) + geom_line(aes(time, fishstock, group = reps), alpha = 0.2) +
 geom_line(aes(time, harvest, group = reps), alpha = 0.05, col="darkgreen")

## Ensemble as a summary ribbon 
stats <- dt[ , mean_sdl(fishstock), by = time]
p1 <- p1 + geom_line(dat=stats, aes(x=time, y=y), col="lightgrey") + 
  geom_ribbon(aes(x = time, ymin = ymin, ymax = ymax),
              fill = "blue", alpha = 0.1, dat=stats)
## Add some reference lines?
p1 <- p1 + geom_abline(intercept=opt$S, slope = 0)          
p1 <- p1 + geom_abline(intercept=xT, slope = 0, lty=2) 

p1

#p2 <- ggplot(dt) + geom_line(aes(time, harvest, group = reps), alpha = 0.2, col="green")
p3 <- ggplot(dt) + geom_line(aes(time, profits, group = reps), alpha = 0.2)
#p2 <- ggplot(dt) + geom_point(aes(time, harvest), alpha=.2)



p0 <- ggplot(subset(dt,reps==sample(1:100))) + 
  geom_line(aes(time, fishstock)) + 
  geom_line(aes(time, harvest), col="darkgreen") + 
  geom_line(aes(time, unharvested), col="lightblue", alpha=.7) +   
  geom_line(aes(time, escapement), col="darkred", alpha=.4) + 
  geom_abline(intercept=opt$S, slope = 0, lwd=1, lty=2)          
p0


## Shows clearly that groups are determined by how many times they got to harvest.  
p4 <- ggplot(dt, aes(total.profit, fill=crashed)) + geom_density(alpha=.8)


## Add discrete classes by total profit
quantile_me <- function(x, ...){
  q <- quantile(x, ...)
  class <- character(length(x))
  for(i in 1:length(q))
    class[x > q[i] ] <- i
  class
}
q <- data.table(reps=total_profit$reps, quantile=quantile_me(total_profit$V1))
setkey(q, reps)
dt <- dt[q]


# fast way to get p0
#simple <- melt(sims, "time")  
#p0 <- ggplot(subset(simple, L1==sample(1:100))) + geom_line(aes(time, value, color=variable))

