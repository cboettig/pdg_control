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

#set.seed(1)
# load my script defining the SDP functions:
# determine_SDP_matrix, find_dp_optim, & ForwardSimulate
source("stochastic_dynamic_programming.R")

# Define all parameters 
delta <- 0.1      # economic discounting rate
OptTime <- 50     # stopping time
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


#' define a profit function, price minus cost
#' @param x is a the grid of state values (profit will evaluate at each of them)
#' @param h is the current harvest level being considered by the algorithm
#' @param p price of fish (Note, optimal will scrap xT if price is high enough!) 
#' @param c fishing extraction costs (per unit effort)
profit <- function(x_grid,h_i, p=1, c=.001){
  ## Havest-based control; havest cannot exceed population size
#  harvest <- sapply(x_grid, function(x_i) min(h_i, x_i))
#  out <- sapply(harvest, function(x) max(0,p*x - c/x))

  ## Effort-based control: 
#  out <- sapply(x_grid, function(x) max(0,p*h_i*x - c*h_i))

  # Another effort-based control cost-function (?)
  harvest <- x_grid*h_i # cpue proportional to population size
  out <- sapply(harvest, function(x) max(0,p*x - c/x))
  out
}

# Set up the grid 
x_grid <- seq(0, 2*K, length=gridsize)  # population size
#h_grid <- x_grid  # vector of havest levels, use same res as stock
h_grid <- seq(0, 2, length=gridsize) # Myers model based on effort!

# Calculate the transition matrix 
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma)

# Find the optimum by dynamic programming 
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta)

## What if we've assumed the wrong model?  
## Estimate a logistic model from the data



## Example plot the results of a single run, against unharvested version  
out <- ForwardSimulate(f, pars, x_grid, h_grid, sigma, x0, opt$D)
dat <- melt(out, id="time")
p0 <- ggplot(dat, aes(time, value, color=variable)) + geom_line() +  
  geom_abline(intercept=opt$S, slope=0, col="black") + # Reed's S,
  geom_abline(intercept = xT, slope=0, col="darkred", lty=3) # unfished Allee 
p0 <- p0 + geom_abline(intercept=e_star, slope=0, col="green", lty = 2) # tippt

#######################################################################
# Now we'll simulate this process many times under this optimal havest#
#######################################################################
sims <- lapply(1:100, function(i){
# simulate the optimal routine on a stoch realization of growth dynamics
    ForwardSimulate(f, pars, x_grid, h_grid, sigma, x0, opt$D)
})
dat <- melt(sims, id="time")

# some stats, geom-smooth should do this?
m <- cast(dat, time ~ variable, mean) # mean population
err <- cast(dat, time ~ variable, sd) # sd population
m_f <- m$fishstock
err_f <- err$fishstock 
m_h <- m$harvest
err_h <- err$harvest 

p1 <- ggplot(dat) +
      # Replicate harvested dynamics
      geom_line(aes(time, value, group = L1), data = 
                subset(dat, variable == "fishstock"), alpha = 0.4) + 
      ## Mean & SD for population
      geom_ribbon(aes(x = time, ymin = m_f - err_f, ymax = m_f + err_f),
                  fill = "darkblue", alpha = 0.7)  +
      geom_line(aes(time, m_f), col = "lightblue")  +
#      geom_abline(intercept=opt$S, slope = 0) + # Reed's S: optimal escapement 
      geom_abline(intercept=xT, slope = 0, lty=2) + # Allee threshold
      ## And the same for harvest levels
      geom_line(aes(time, value, group = L1), 
            data = subset(dat, variable == "harvest"),  alpha=.2, col = "darkgreen") +
      geom_ribbon(aes(x = time, ymin = m_h - err_h, ymax = m_h + err_h),
                  fill = "darkgreen", alpha = 0.4)  +
      geom_line(aes(time, m_h), col = "lightgreen")  +
      geom_abline(intercept = e_star, slope = 0, col = "lightgreen", lwd=1,lty=2) 

optimal_crashed = subset(dat, variable == "fishstock" & 
                         time == OptTime-1 & value < xT)
p1 <- p1 + opts(title=sprintf("Optimal Harvest dynamics, %d populations crash",
  dim(optimal_crashed)[1]))
print(p1)

## make the crashed trajectories stand out?
#p1 <- p1 + geom_line(aes(time, value, group = L1), 
#          data = subset(dat, variable == "harvest" & 
#          (L1 %in% optimal_crashed$L1)),  col = "darkgreen", alpha = 0.5) +
#     geom_line(aes(time, value, group = L1), 
#          data = subset(dat, variable == "fishstock" &
#          (L1 %in% optimal_crashed$L1)),  col = "darkblue", alpha = 0.5) 


crashed = subset(dat, variable =="unharvested" & time == OptTime-1 & value < xT)
p2 <- ggplot(dat) + 
  geom_line(aes(time, value, group = L1), 
            data = subset(dat, variable == "unharvested"), alpha=.2) + 
  geom_line(aes(time, cast(dat, time ~ variable, mean)$unharvested)) 

p2 <- p2 + opts(title=sprintf("Unfished dynamics, %d populations crash",
  dim(crashed)[1]))



#ggsave("samplerun.png", plot=p0)
#ggsave("fished.png", plot=p1)
#ggsave("unfished.png", plot=p2)








# check out the deterministic dynamics
x <- numeric(OptTime)
x[1] <- K
for(t in 1:(OptTime-1))
  x[t+1] <- f(x[t],e_star+.01,pars)

