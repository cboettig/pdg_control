




# Reed Model
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

 Implements a numerical version of the SDP described in:
 
   Sethi, G., Costello, C., Fisher, A., Hanemann, M., and Karp, L. (2005). 
   Fishery management under multiple uncertainty. Journal of Environmental
   Economics and Management, 50(2), 300-318. doi:10.1016/j.jeem.2004.11.005

   Reed, W.J., 1979. Optimal Escapement Levels in Stochastic
   and Deterministic Harvesting Models. Journal of Environmental 
   Economics and Management. 6: 350-363.


Clear the workspace and load package dependencies: 




### Define all parameters 


```r
delta <- 0.1      # economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
reward <- 1       # bonus for satisfying the boundary condition
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
control = "harvest"         # control variable is total harvest, h = e * x
price <- 1
cost <- .12
```




And a state equation with an allee effect


```r
f <- Myer_harvest
pars <- c(1, 2, 6) 
p <- pars # shorthand 
K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
xT <- p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 # allee threshold
e_star <- (p[1] * sqrt(p[3]) - 2) / 2 ## Bifurcation point 
control <- "harvest"          # control variable can be harvest or effort 
price <- 1
cost <- .01
```






Our initial condition is the equilibrium size (note the stochastic deflation of mean)


```r
x0 <- K - sigma_g ^ 2 / 2 
```




and we use a harvest-based profit function with default parameters


```r
profit <- profit_harvest(p=price, c = cost) 
```





Set up the discrete grids for stock size and havest levels (which will use same resolution as for stock). 



```r
x_grid <- seq(0, 2 * K, length = gridsize)  
h_grid <- x_grid  
```





### Calculate the transition matrix (with noise in growth only)      
We calculate the stochastic transition matrix for the probability of going from any state \(x_t \) to any other state \(x_{t+1}\) the following year, for each possible choice of harvest \( h_t \).  This provides a look-up table for the dynamic programming calculations. 



```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
```







```r
# calculate the transition matrix by simulation, generic to types of noise
require(snowfall) 
sfInit(parallel=TRUE, cpu=4)
SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps=999)
```




### Find the optimum by dynamic programming 
Bellman's algorithm to compute the optimal solution for all possible trajectories.


```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
```




### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above.


```r
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})
```







## Summarize and plot the results                                                   
Make data tidy (melt), fast (data.tables), and nicely labeled.


```r
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
```




### Plots 

Let's begin by looking at the dynamics of a single replicate. The line shows Reed's S, the level above which the stock should be harvested (where catch should be the difference between stock and S).  To confirm that this policy is being followed, note that harvesting only occurs when the stock is above this line, and harvest is proportional to the amount by which it is above. 


```r
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, fishstock)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_line(aes(time, harvest), col="darkgreen") 
```

![plot of chunk unnamed-chunk-11](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-ex-out-unnamed-chunk-112.png) 



This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again, along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). 


```r
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) + 
  geom_abline(intercept=xT, slope = 0, lty=2) 
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
```

![plot of chunk unnamed-chunk-12](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-ex-out-unnamed-chunk-121.png) 


We can also look at the harvest dynamics:


```r
p1 + geom_line(aes(time, harvest, group = reps), alpha = 0.1, col="darkgreen")
```

![plot of chunk unnamed-chunk-13](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-ex-out-unnamed-chunk-131.png) 


This strategy is supposed to be a constant-escapement strategy. We can visualize the escapement: 


```r
p1 + geom_line(aes(time, escapement, group = reps), alpha = 0.1, col="darkgrey")
```

![plot of chunk unnamed-chunk-14](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-ex-out-unnamed-chunk-141.png) 




### Computing additional statistics about the data
In this section we add some additional information to our data.table on the profits obtained by each replicate.  The algorithm has supposedly maximized the expected profit, so it is useful to look at both the mean total profit and the distribution.  Despite this maximization, the distribution can be rather lop-sided or even bimodal. 

Which replicates crashed?  Which met the boundary requirment and recieved the reward value at the end?


```r
crashed <- dt[time==OptTime, fishstock == 0, by=reps]
rewarded <- dt[time==OptTime, fishstock > xT, by=reps]
```




Let's compute the profits at each time-step for each replicate. 
Using `data.table` to evaluate our profit function over the stock and harvest levels requires indexing our data:



```r
dt <- data.table(dt, id=1:dim(dt)[1])
profits <- dt[, profit(fishstock, harvest), by=id]
```




Merging this calculation back into our data table using fast join (needs to define 'id' as a key on which to match things up though). 


```r
setkey(dt, id)
setkey(profits, id)
dt <- dt[profits]
setnames(dt, "V1", "profits")
setkey(dt, reps)
```




Compute total profit by summing over each timeseries (including the reward for satisfying the terminal boundary condition, if any). 



```r
total_profit <- dt[,sum(profits), by=reps]
total_profit <- total_profit + rewarded$V1 * reward 
```





Add these three columns to the data.table (fast join and re-label):


```r
setkey(total_profit, reps)
setkey(crashed, reps)
setkey(rewarded, reps)
dt <- dt[total_profit]
dt <- dt[crashed]
dt <- dt[rewarded]
setnames(dt, c("V1", "V1.1", "V1.2"), c("total.profit", "crashed", "rewarded"))
```






#### Profit plots
Since the optimal strategy maximizes expected profit, it may be more useful to look at the distribution statistics of profit over time:


```r
stats <- dt[ , mean_sdl(profits), by = time]
p1 + geom_line(dat=stats, aes(x=time, y=y), col="lightgrey") + 
  geom_ribbon(aes(x = time, ymin = ymin, ymax = ymax),
              fill = "darkred", alpha = 0.2, dat=stats)
```

![plot of chunk unnamed-chunk-20](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-ex-out-unnamed-chunk-20.png) 



Total profits


```r
ggplot(dt, aes(total.profit, fill=crashed)) + geom_histogram(alpha=.8)
```

![plot of chunk unnamed-chunk-21](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-ex-out-unnamed-chunk-21.png) 



#### Add discrete classes by total profit

Sometimes I'd like to color code the replicates by profit, to see if there are particular patterns in stock dynamics of the most profitable and least profitable lines.  Adding discrete profit classes to the data table makes this possible:


```r
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
```




Then we can plot the fishstock trajectories, indicating which derive the highest and smallest profits by color code: 


```r
ggplot(subset(dt, quantile %in% c(1,4))) + 
  geom_line(aes(time, fishstock, group = reps, color=quantile), alpha = 0.6) 
```

![plot of chunk unnamed-chunk-23](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-ex-out-unnamed-chunk-23.png) 



### Visualizing the optimal policy
Note that when the boundary is sufficiently far away, i.e. for the first couple timesteps, the optimal policy is stationary.  The optimal policy is shown here over time, where the color indicates the harvest recommended for each possible stock value at that time (shown on the vertical axis).  Note that below a certain stock value, harvesting is not recommended and the dots turn red (Reed's constant escapement rule!)  However, at very low values, harvesting starts again (orange dots), because of the allee effect - these populations are doomed anyway, so may as well fish all that remains.

Note that interestingly, populations just below the allee threshold are given the chance to be rescued stochastically early on - that small chance that they recover is worth the expected loss.  The "no-harvest" zones stand out clearly in the red areas of this graph.


```r
policy <- melt(opt$D)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
p5 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p5
```

![plot of chunk policy](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-ex-out-policy.png) 


The harvest intensity is limited by the stock size.  If instead we look at the difference between proposed harvest intensity and stock,
then the red zones correspond to places where harvest equals stock, i.e. we slam the population. The allee band is clearly seen. Just for fun, we overlay the replicate dynamics on this plot 



```r
p6 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p6 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1, data=dt)
```

![plot of chunk policy2](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-ex-out-policy2.png) 

