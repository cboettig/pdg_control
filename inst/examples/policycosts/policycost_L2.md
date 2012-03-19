




# L2 Policy Costs 
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

## Setup the system



```r
rm(list=ls())   
require(pdgControl)
require(reshape2)
require(ggplot2)
require(data.table)
```



This example illustrates the impact of adding a cost to changing the harvest level between years 

### Define all parameters 


```r
delta <- 0.01     # SMALLER economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
reward <- 0       # bonus for satisfying the boundary condition
```




we'll use log normal noise functions


```r
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
```





Chose the state equation / population dynamics function


```r
f <- BevHolt                # Select the state equation
pars <- c(2, 4)             # parameters for the state equation
K <- (pars[1] - 1)/pars[2]  # Carrying capacity 
xT <- 0                     # boundary conditions
```




Our initial condition is the equilibrium size (note the stochastic deflation of mean)


```r
x0 <- K - sigma_g ^ 2 / 2 
```




and we use a harvest-based profit function with default parameters


```r
profit <- profit_harvest(price_fish = 1, cost_stock_effect = 0,
 operating_cost = 0.1)
```




Set up the discrete grids for stock size and havest levels


```r
x_grid <- seq(0, 1.2 * K, length = gridsize)  
h_grid <- seq(0, 0.8 * K, length = gridsize)  
```




### Calculate the stochastic transition matrix
We calculate the stochastic transition matrix for the probability of going from any state \(x_t \) to any other state \(x_{t+1}\) the following year, for each possible choice of harvest \( h_t \).  This provides a look-up table for the dynamic programming calculations. Note that this only includes uncertainty in the growth rate (projected stock next year). 


```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
```



### Find the optimum by dynamic programming 
We use Bellman's dynamic programming algorithm to compute the optimal solution for all possible trajectories, ignoring potential policy costs as before.  We will later use this solution to compare against the optimal solution with policy costs.


```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
```




A modified algorithm lets us include a penalty of magnitude `P` and a functional form that can be an `L1` norm, `L2`  norm, `asymmetric` L1 norm (costly to lower harvest rates), fixed cost, or `none` (no cost).  Here is an asymmetric norm example.  Note that this calculation is considerably slower. 


```r
policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                    profit, delta, reward, P = 5, penalty = "L2")
```





### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above.  We use a modified simulation function that can simulate an alternate policy (the Reed optimum, where policy costs are zero, `opt$D` ) and a focal policy, `policycost$D`



```r
sims <- lapply(1:100, function(i)
  simulate_optim(f, pars, x_grid, h_grid, x0, policycost$D, z_g, z_m, z_i, opt$D)
  )
```





Make data tidy (melt), fast (data.tables), and nicely labeled.


```r
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
```




### Plots 

A single replicate, alternate dynamics should show the Reed optimum, while harvest/fishstock should show the impact of having policy costs. 


```r
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, alternate)) +
  geom_line(aes(time, fishstock), col="darkblue") +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_line(aes(time, harvest), col="purple") + 
  geom_line(aes(time, harvest_alt), col="darkgreen") 
```

![plot of chunk rep1](http://farm8.staticflickr.com/7038/6982916703_4b9b8aed90_o.png) 



We can visualize the equilibrium policy for each possible harvest:



```r
policy <- sapply(1:length(h_grid), function(i) policycost$D[[i]][,1])
ggplot(melt(policy)) + 
  geom_point(aes(h_grid[Var2], (x_grid[Var1]), col=h_grid[value] - h_grid[Var2])) + 
    labs(x = "prev harvest", y = "fishstock") +
      scale_colour_gradientn(colours = rainbow(4)) 
```

![plot of chunk unnamed-chunk-2](http://farm8.staticflickr.com/7179/6982916875_8f8b28c3d6_o.png) 


Here we plot previous harvest against the recommended harvest, coloring by stocksize.  Note this swaps the y axis from above with the color density.  Hence each x-axis value has all possible colors, but they map down onto a subset of optimal harvest values (depending on their stock). 


```r
policy <- sapply(1:length(h_grid), function(i) policycost$D[[i]][,1])
ggplot(melt(policy)) + 
  geom_point(aes(h_grid[Var2], (h_grid[value]), col = x_grid[Var1]), position=position_jitter(w=.005,h=.005), alpha=.5) + 
    labs(x = "prev harvest", y = "harvest") +
      scale_colour_gradientn(colours = rainbow(4)) 
```

![plot of chunk unnamed-chunk-3](http://farm8.staticflickr.com/7044/6982917077_abecdd2da6_o.png) 



Against the policy with no cost (shown over time) 


```r
policy <- melt(opt$D)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$alternate) )
ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col= h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) + 
  geom_abline(intercept=opt$S, slope = 0) 
```

![plot of chunk no_policy_cost_vis](http://farm8.staticflickr.com/7069/6982917243_cdf6b5965b_o.png) 


### Profits



```r
dt <- data.table(dt, id=1:dim(dt)[1])
profits <- dt[, profit(fishstock, harvest), by=id]
```




Merge in profits to data.table (should be a way to avoid having to do these joins?)


```r
setkey(dt, id)
setkey(profits, id)
dt <- dt[profits]
setnames(dt, "V1", "profits")
```




merge in total profits to data.table


```r
total_profit <- dt[,sum(profits), by=reps]
setkey(total_profit, reps)
setkey(dt, reps)
dt <- dt[total_profit]
setnames(dt, "V1", "total.profit")
```






```r
ggplot(dt, aes(total.profit)) + geom_histogram(alpha=.8)
```

![plot of chunk unnamed-chunk-7](http://farm8.staticflickr.com/7058/6836790866_d942e6e3a1_o.png) 




```r
save(list=ls(), file="L2.rda")
```




The mean dynamics of the state


```r
stats <- dt[ , mean_sdl(fishstock), by = time]
ggplot(stats) +   geom_ribbon(aes(x = time, ymin = ymin, ymax = ymax),
                fill = "darkblue", alpha = 0.2, dat=stats) +
                geom_line(aes(x=time, y=y), lwd=1) 
```

![plot of chunk unnamed-chunk-9](http://farm8.staticflickr.com/7177/6982917599_df2afaac5d_o.png) 


The mean dynamics of the control


```r
stats <- dt[ , mean_sdl(harvest), by = time]
ggplot(stats) +   geom_ribbon(aes(x = time, ymin = ymin, ymax = ymax),
                fill = "darkblue", alpha = 0.2, dat=stats) +
                geom_line(aes(x=time, y=y), lwd=1) 
```

![plot of chunk unnamed-chunk-10](http://farm8.staticflickr.com/7058/6836791324_b5e91c2053_o.png) 


