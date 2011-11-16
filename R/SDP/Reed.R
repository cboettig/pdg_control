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
# Dependencies:
#   functions in scripts: 
#     stochastic_dynamic_programming.R, 
#     population_models.R
#   Requires library: "ggplot2" (plotting)
#   Recommends library: "snowfall" (for simulated transition matrix)
#
#
# Fish population dynamics:
# X_{t+1} = Z_n f(X_n) 


rm(list=ls())   # Start wtih clean workspace 
require(ggplot2) # nicer plotting package

# load my script defining the SDP functions:
# determine_SDP_matrix, find_dp_optim, & ForwardSimulate
source("stochastic_dynamic_programming.R")

# load the library of population models
source("population_models.R")

# Define all parameters 
delta <- 0.1      # economic discounting rate
OptTime <- 50     # stopping time
sigma_g <- 0.2    # Noise in population growth
gridsize <- 100   # gridsize (discretized population)
sigma_m <- .0     # noise in stock assessment measurement
sigma_i <- .0     # noise in implementation of the quota
interval <- 1     # period of updating the stock assessment



# Chose the state equation / population dynamics function
f <- BevHolt              # Select the state equation
pars <- c(2,4)            # parameters for the state equation
K <- (pars[1]-1)/pars[2]  # Carrying capacity 
xT <- 0                   # boundary conditions
e_star <- 0               # model's bifurcation point (just for reference)
 control = "harvest"         # control variable is total harvest, h = e * x

## An alternative state equation, with allee effect: (uncomment to select)
#f <- Myer
#pars <- c(1, 2, 6) 
#p <- pars # shorthand 
#K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
#xT <- p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 # allee threshold
#e_star <- (p[1] * sqrt(p[3]) - 2) / 2 ## Bifurcation point 
#control <- "effort"          # control variable is harvest effort, e = h / x

x0 <- K # initial condition

#' Define a profit function, price minus cost
#' @param x is a the grid of state values (profit will evaluate at each of them)
#' @param h is the current harvest level being considered by the algorithm
#' @param p price of fish (Note, optimal will scrap xT if price is high enough!) 
#' @param c fishing extraction costs (per unit effort)
profit <- function(x_grid, h_i, p = 1, c = 0.001, type=control){
  if(type=="harvest")
    harvest <- sapply(x_grid, function(x_i) min(h_i, x_i)) # Harvest-based control
  else if(type=="effort")
    harvest <- x_grid * h_i # Effort-based control 
  sapply(harvest, function(x) max(0, p * x - c / x))
}

# Set up the discrete grids
x_grid <- seq(0, 2 * K, length = gridsize)  # population size
#h_grid <- x_grid  # vector of havest levels, use same resolution as for stock
h_grid <- seq(0, 2, length=gridsize) 


#######################################################################
# Calculate the transition matrix (with noise in growth only)         #
#######################################################################
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g)

## calculate the transition matrix by simulation 
#require(snowfall) # use parallelization since this can be slow
#sfInit(parallel=TRUE, cpu=4)
#SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, sigma_g, sigma_m, sigma_i, reps=99)


#######################################################################
# Find the optimum by dynamic programming                             #
#######################################################################
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, reward=100)


# What if parameter estimation is inaccurate? e.g.:
#pars[1] <- pars[1] * 0.95


#######################################################################
# Now we'll simulate this process many times under this optimal havest#
#######################################################################
sims <- lapply(1:100, function(i){
# simulate the optimal routine on a stoch realization of growth dynamics
  ForwardSimulate(f, pars, x_grid, h_grid, sigma_g, x0, opt$D, 
                    sigma_m, sigma_i, interval=1)
})




########################################################################
#  Summarize and plot data                                             #
########################################################################
dat <- melt(sims, id="time") # reshapes the data matrix to "long" form

# some stats on the replicates, ((geom_smooth should do this?))
m <- cast(dat, time ~ variable, mean) # mean population
err <- cast(dat, time ~ variable, sd) # sd population
m_f <- m$fishstock
err_f <- err$fishstock 
m_h <- m$harvest
err_h <- err$harvest 


## Creat the plot
p1 <- ggplot(dat) +
      # Replicate harvested dynamics
      geom_line(aes(time, value, group = L1), data = 
                subset(dat, variable == "fishstock"), alpha = 0.2) + 
      ## Mean & SD for population
      geom_ribbon(aes(x = time, ymin = m_f - err_f, ymax = m_f + err_f),
                  fill = "darkblue", alpha = 0.4)  +
      geom_line(aes(time, m_f), col = "lightblue")  +
#      geom_abline(intercept=opt$S, slope = 0) + # show Reed's S: optimal escapement 
      geom_abline(intercept=xT, slope = 0, lty=2) + # show Allee threshold
      ## And the same for harvest-effort levels
#      geom_line(aes(time, value, group = L1), 
#            data = subset(dat, variable == "harvest"),  alpha=.2, col = "darkgreen") +
      geom_ribbon(aes(x = time, ymin = m_h - err_h, ymax = m_h + err_h),
                  fill = "darkgreen", alpha = 0.4)  +
      geom_line(aes(time, m_h), col = "lightgreen")  #+
#     geom_abline(intercept = e_star, slope = 0, col = "lightgreen", lwd=1,lty=2) 


## Count how many crashed and add it in a plot title
optimal_crashed = subset(dat, variable == "fishstock" & 
                         time == OptTime-1 & value < xT)
p1 <- p1 + opts(title = sprintf("Optimal Harvest dynamics, %d populations crash",
                                dim(optimal_crashed)[1]))
print(p1)

## extra plots are avialable in plots.R, inculding unharvested dynamics,
## plot of a single harvested replicate, and plot of the profit over time.  
# source("plots.R")
