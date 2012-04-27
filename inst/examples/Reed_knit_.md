<!--roptions dev='png', fig.width=7, fig.height=5, fig.path='figure/', tidy=FALSE, warning=FALSE, comment=NA, message=FALSE, refresh=2-->

<!--begin.rcode echo=FALSE 
render_markdown()
opts_knit$set(upload = TRUE)
require(socialR)
options(flickrOptions=list(
  description="https://github.com/cboettig/pdg_control/blob/master/inst/examples/",
  tags="stochpop, pdg_control"))
opts_knit$set(upload.fun = flickr.url)
options(device = function(width = 5, height = 5) {
    pdf(NULL, width = width, height = height)
})
end.rcode-->

# Reed Model
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

 Implements a numerical version of the SDP described in:

   Reed, W.J., 1979. Optimal Escapement Levels in Stochastic
   and Deterministic Harvesting Models. Journal of Environmental 
   Economics and Management. 6: 350-363.


Clear the workspace and load package dependencies: 
<!--begin.rcode setup, echo=FALSE
rm(list=ls())   
require(pdgControl)
require(reshape2)
require(ggplot2)
require(data.table)
end.rcode-->


### Define all parameters 
<!--begin.rcode
delta <- 0.1      # economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
reward <- 1       # bonus for satisfying the boundary condition
end.rcode-->

we'll use log normal noise functions. 
For Reed, only `z_g` will be random.
Sethi et al will add the other forms
<!--begin.rcode noise_dists
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
end.rcode-->


Chose the state equation / population dynamics function
And a state equation with an allee effect
<!--begin.rcode
f <- BevHolt
pars <- c(1.5, 5)
K <- pars[2]/(1-pars[1])
xT <- 0
price <- 1
cost <- .01
end.rcode-->



Our initial condition is the equilibrium size (note the stochastic deflation of mean)
<!--begin.rcode Xo
x0 <- K - sigma_g ^ 2 / 2 
end.rcode-->

and we use a harvest-based profit function with default parameters
<!--begin.rcode profit
profit <- profit_harvest(price=price, c0 = cost) 
end.rcode-->


Set up the discrete grids for stock size and havest levels (which will use same resolution as for stock). 

<!--begin.rcode grid
x_grid <- seq(0, 2 * K, length = gridsize)  
h_grid <- x_grid  
end.rcode-->


### Calculate the transition matrix (with noise in growth only)      
We calculate the stochastic transition matrix for the probability of going from any state \(x_t \) to any other state \(x_{t+1}\) the following year, for each possible choice of harvest \( h_t \).  This provides a look-up table for the dynamic programming calculations. 

<!--begin.rcode SDP_Mat
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
end.rcode-->


<!--begin.rcode parallel, eval=FALSE, include=FALSE
# calculate the transition matrix by simulation, generic to types of noise
require(snowfall) 
sfInit(parallel=TRUE, cpu=4)
SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps=999)
end.rcode-->

### Find the optimum by dynamic programming 
Bellman's algorithm to compute the optimal solution for all possible trajectories.
<!--begin.rcode find_dp_opt 
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
end.rcode-->

### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above.
<!--begin.rcode simulate 
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})
end.rcode-->




## Summarize and plot the results                                                   
Make data tidy (melt), fast (data.tables), and nicely labeled.
<!--begin.rcode tidy
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
end.rcode-->

### Plots 

Let's begin by looking at the dynamics of a single replicate. The line shows Reed's S, the level above which the stock should be harvested (where catch should be the difference between stock and S).  To confirm that this policy is being followed, note that harvesting only occurs when the stock is above this line, and harvest is proportional to the amount by which it is above. 

<!--begin.rcode p0
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, fishstock)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_line(aes(time, harvest), col="darkgreen") 
end.rcode-->


This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again, along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). 

<!--begin.rcode p1
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) + 
  geom_abline(intercept=xT, slope = 0, lty=2) 
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
end.rcode-->

We can also look at the harvest dynamics:
<!--begin.rcode p2
p1 + geom_line(aes(time, harvest, group = reps), alpha = 0.1, col="darkgreen")
end.rcode-->

This strategy is supposed to be a constant-escapement strategy. We can visualize the escapement: 

<!--begin.rcode p3
p1 + geom_line(aes(time, escapement, group = reps), alpha = 0.1, col="darkgrey")
end.rcode-->



### Computing additional statistics about the data
In this section we add some additional information to our data.table on the profits obtained by each replicate.  The algorithm has supposedly maximized the expected profit, so it is useful to look at both the mean total profit and the distribution.  Despite this maximization, the distribution can be rather lop-sided or even bimodal. 

Which replicates crashed?  Which met the boundary requirment and recieved the reward value at the end?
<!--begin.rcode crashed
crashed <- dt[time==OptTime, fishstock == 0, by=reps]
rewarded <- dt[time==OptTime, fishstock > xT, by=reps]
end.rcode-->


Add these three columns to the data.table (fast join and re-label):
<!--begin.rcode
setkey(crashed, reps)
setkey(rewarded, reps)
dt <- dt[crashed]
dt <- dt[rewarded]
setnames(dt, c("V1.1", "V1.2"), c("crashed", "rewarded"))
end.rcode-->



#### Profit plots
Since the optimal strategy maximizes expected profit, it may be more useful to look at the distribution statistics of profit over time:
<!--begin.rcode profitplot
stats <- dt[ , mean_sdl(profits), by = time]
p1 + geom_line(dat=stats, aes(x=time, y=y), col="lightgrey") + 
  geom_ribbon(aes(x = time, ymin = ymin, ymax = ymax),
              fill = "darkred", alpha = 0.2, dat=stats)
end.rcode-->


Total profits
<!--begin.rcode totalprofit
ggplot(dt, aes(total.profit, fill=crashed)) + geom_histogram(alpha=.8)
end.rcode-->


#### Add discrete classes by total profit

Sometimes I'd like to color code the replicates by profit, to see if there are particular patterns in stock dynamics of the most profitable and least profitable lines.  Adding discrete profit classes to the data table makes this possible:
<!--begin.rcode quantiles
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
end.rcode-->

Then we can plot the fishstock trajectories, indicating which derive the highest and smallest profits by color code: 
<!--begin.rcode colorquantiles
ggplot(subset(dt, quantile %in% c(1,4))) + 
  geom_line(aes(time, fishstock, group = reps, color=quantile), alpha = 0.6) 
end.rcode-->


### Visualizing the optimal policy
Note that when the boundary is sufficiently far away, i.e. for the first couple timesteps, the optimal policy is stationary.  The optimal policy is shown here over time, where the color indicates the harvest recommended for each possible stock value at that time (shown on the vertical axis).  Note that below a certain stock value, harvesting is not recommended and the dots turn red (Reed's constant escapement rule!)  However, at very low values, harvesting starts again (orange dots), because of the allee effect - these populations are doomed anyway, so may as well fish all that remains.

Note that interestingly, populations just below the allee threshold are given the chance to be rescued stochastically early on - that small chance that they recover is worth the expected loss.  The "no-harvest" zones stand out clearly in the red areas of this graph.

<!--begin.rcode policy, fig.width=9
policy <- melt(opt$D)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
p5 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p5
end.rcode-->

The harvest intensity is limited by the stock size.


<!--begin.rcode policy2, fig.width=10
p6 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p6 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1, data=dt)
end.rcode-->

