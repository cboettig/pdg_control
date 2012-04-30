`ro cache=TRUE, tidy=FALSE, warning=FALSE, comment=NA, message=FALSE, refresh=1 or`

``` {r  echo=FALSE }
opts_knit$set(upload.fun = socialR::flickr.url)
options(device = function(width = 5, height = 5) {
    pdf(NULL, width = width, height = height)
})
````

# Reed Model
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

 Implements a numerical version of the SDP described in:

   Reed, W.J., 1979. Optimal Escapement Levels in Stochastic
   and Deterministic Harvesting Models. Journal of Environmental 
   Economics and Management. 6: 350-363.


Clear the workspace and load package dependencies: 
``` {r  setup, echo=FALSE}
rm(list=ls())   
require(pdgControl)
require(reshape2)
require(ggplot2)
require(data.table)
````


### Define parameters 

``` {r }
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
````

we'll use log normal noise functions. 
For Reed, only `z_g` will be random.
Sethi et al will add the other forms

``` {r  noise_dists}
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
````


Chose the state equation / population dynamics function

``` {r }
f <- BevHolt
````

Note that the `pdg_control` pacakge already has a definition for the `BevHolt` function, (typing the function name prints the function)

``` {r }
BevHolt
````

That is, \( f(x,h) = \frac{A x}{1 + B x} \)

Of course we could pass in any custom function of stocksize `x`, harvest `h` and parameter vector `p` in place of `BevHolt`.  Note that we would need to write this function explicitly so that it can take vector values of `x` (i.e. uses `sapply`), an annoying feature of `R` for users comming from Matlab.  


We must now define parameters for the function.  Note that the positive stationary root of the model is given by \( (A-1)/B \), which we'll store for future reference as `K`.  

``` {r }
pars <- c(1.5, 0.05)
K <- (pars[1] - 1)/pars[2]
````



and we use a harvest-based profit function with default parameters

``` {r  profit}
profit <- profit_harvest(price=1, c0 = 0.01) 
````

The profit_harvest function has the form \( \Pi = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x} \), conditioned on \( h > x \) and \(x > 0 \).  Note that the R code defines a function from another function using a trick known as a _closure_.  Again we could write a custom profit function as long as it can take a vector stock size `x` and a scalar harvest level `h`.  Details for provided functions can be found in the manual, i.e. `?profit_harvest`. 


Now we must set up the discrete grids for stock size and havest levels (which will use same resolution as for stock), in order to calculate the SDP solution.   Here we set the gridsize to 100.  

``` {r  grid}
x_grid <- seq(0, 2 * K, length = 100)  
h_grid <- x_grid  
````


### Calculate the transition matrix (with noise in growth only)      
We calculate the stochastic transition matrix for the probability of going from any state \(x_t \) to any other state \(x_{t+1}\) the following year, for each possible choice of harvest \( h_t \).  This provides a look-up table for the dynamic programming calculations.  Take a look at the R code and the documentation for `determine_SDP_matrix` to see what this function is actually doing.  The implementation is quite simple.   

``` {r  SDP_Mat}
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
````


### Find the optimum by dynamic programming

Bellman's algorithm to compute the optimal solution for all possible trajectories.

``` {r  find_dp_opt }
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
````

Note that `SDP_Mat` is specified from the calculation above, as are our grids and our profit function. `OptTime` is the stopping time.  `xT` specifies a boundary condition at the stopping time. A reward for meeting this boundary must be specified for it to make any difference.  `delta` indicates the economic discount rate. Again, details are in the function documentation.   


In the Sethi case, computing the distribution over multiple sources of noise is actually quite difficult.  Simulation turns out to be more efficient than numerically integrating over each distribution.  This code parallelizes the operation over four cores, but can be scaled to an arbitrary cluster. Since we're focused on the Reed example for the moment, we can ignore this step.   

``` {r  parallel, eval=FALSE, include=FALSE}
require(snowfall) 
sfInit(parallel=TRUE, cpu=4)
SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps=999)
````

### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the Reed optimal harvest policy determined above.

``` {r  simulate }
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i)
})
````

The forward simulation algorithm needs an initial condition `x0` which we set equal to the carrying capacity, as well as our population dynamics `f`, parameters `pars`, grids, and noise coefficients.  Recall in the Reed case only `z_g`, growth, is stochastic.  


## Summarize and plot the results                                                   

R makes it easy to work with this big replicate data set.  We make data tidy (melt), fast (data.tables), and nicely labeled.

``` {r  tidy}
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
````

### Plots 

Let's begin by looking at the dynamics of a single replicate. The line shows Reed's S, the level above which the stock should be harvested (where catch should be the difference between stock and S).  To confirm that this policy is being followed, note that harvesting only occurs when the stock is above this line, and harvest is proportional to the amount by which it is above.  Change the replicate `reps==` to see the results from a different replicate.  

``` {r  p0}
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, fishstock)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_line(aes(time, harvest), col="darkgreen") 
````


This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again, along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). 

``` {r  p1}
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) + 
  geom_abline(intercept=xT, slope = 0, lty=2) 
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
````

We can also look at the harvest dynamics:

``` {r  p2}
p1 + geom_line(aes(time, harvest, group = reps), alpha = 0.1, col="darkgreen")
````

This strategy is supposed to be a constant-escapement strategy. We can visualize the escapement: 

``` {r  p3}
p1 + geom_line(aes(time, escapement, group = reps), alpha = 0.1, col="darkgrey")
````





### Visualizing the optimal policy
Note that when the boundary is sufficiently far away, i.e. for the first couple timesteps, the optimal policy is stationary.  The optimal policy is shown here over time, where the color indicates the harvest recommended for each possible stock value at that time (shown on the vertical axis).  Note that below a certain stock value, harvesting is not recommended and the dots turn red (Reed's constant escapement rule!)  However, at very low values, harvesting starts again (orange dots), because of the allee effect - these populations are doomed anyway, so may as well fish all that remains.

Note that interestingly, populations just below the allee threshold are given the chance to be rescued stochastically early on - that small chance that they recover is worth the expected loss.  The "no-harvest" zones stand out clearly in the red areas of this graph.

``` {r  policy, fig.width=9}
policy <- melt(opt$D)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
p5 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p5
````

The harvest intensity is limited by the stock size.


``` {r  policy2, fig.width=10}
p6 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p6 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1, data=dt)
````

