




# Perfect Policy, Imperfect Implementation


## Model setup 

Clear the workspace and load package dependencies: 



Define parameters


```r
delta <- 0.1      # economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
reward <- 1       # bonus for satisfying the boundary condition
```




Use log-normal noise functions


```r
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
```




Chose the state equation / population dynamics function


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




Set up discrete grids for stock size and havest levels (which will use same resolution as for stock). 


```r
x_grid <- seq(0, 2 * K, length = gridsize)  
h_grid <- x_grid  
```





### Calculate the transition matrix (with noise in growth only)      
We calculate the stochastic transition matrix for the probability of going from any state \(x_t \) to any other state \(x_{t+1}\) the following year, for each possible choice of harvest \( h_t \).  This provides a look-up table for the dynamic programming calculations. 


```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
```




### Find the optimum by dynamic programming 
Bellman's algorithm to compute the optimal solution for all possible trajectories.


```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
```




### The optimal policy is implemented imperfectly
We add implementation noise: an imperfect implementation (though a symmetric one -- on average the implementation is not worse than assumed by the optimal solution, it is simply variable). 


```r
sigma_i <- 0.4
```




### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above, but with this additional implementation error


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
This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again, along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). 


```r
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) + 
  geom_abline(intercept=xT, slope = 0, lty=2) 
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
```

![plot of chunk fishstock](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-fishstock6.png) 


We can also look at the harvest dynamics:


```r
p1 + geom_line(aes(time, harvest, group = reps), alpha = 0.1, col="darkgreen")
```

![plot of chunk harvest](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-harvest6.png) 


This strategy is supposed to be a constant-escapement strategy. We can visualize the escapement and see if it is less variable than fish stock, and if it is near Reed's S: 


```r
p1 + geom_line(aes(time, escapement, group = reps), alpha = 0.1, col="darkgrey")
```

![plot of chunk escapement](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-escapement6.png) 


### Computing additional statistics about the data
In this section we add some additional information to our data.table on the profits obtained by each replicate.  The algorithm has supposedly maximized the expected profit, so it is useful to look at both the mean total profit and the distribution.  Despite this maximization, the distribution can be rather lop-sided or even bimodal. 

Which replicates crashed?  Which met the boundary requirment and recieved the reward value at the end?


```r
crashed <- dt[time==OptTime, fishstock == 0, by=reps]
rewarded <- dt[time==OptTime, fishstock > xT, by=reps]
```




A total of \Sexpr{sum(crashed$V1)} crash.

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
```



```
Error: columns of j don't evaluate to consistent types for each group: result for group 51 has column 1 type 'logical' but expecting type 'numeric'
```



```r
p1 + geom_line(dat=stats, aes(x=time, y=y), col="lightgrey") + 
  geom_ribbon(aes(x = time, ymin = ymin, ymax = ymax),
              fill = "darkred", alpha = 0.2, dat=stats)
```



```
Error: object 'stats' not found
```





Total profits


```r
ggplot(dt, aes(total.profit, fill=crashed)) + geom_histogram(alpha=.8)
```

![plot of chunk totals](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-totals6.png) 


## Compare to a non-optimal solution
Compare another model, that likewise assumes no implementation error, and also makes a mistake in its estimate of the growth parameter, making it conservative rather than optimal.



```r
sigma_i <- 0.0
sigma_g <- 0.6
```







```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
```




For the simulated implementation, we add the same implementation error back, and we restore biological growth noise to it's true value


```r
sigma_i <- 0.4
sigma_g <- 0.2
```





### Simulate 
Now we simulate as before


```r
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})
```




## Summarize and plot the results                                                  
Using the code above, recreate the plots for this policy and simulation: 


```r
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
```




### Plots 


```r
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) + 
  geom_abline(intercept=xT, slope = 0, lty=2) 
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
```

![plot of chunk unnamed-chunk-1](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-unnamed-chunk-13.png) 




```r
p1 + geom_line(aes(time, harvest, group = reps), alpha = 0.1, col="darkgreen")
```

![plot of chunk unnamed-chunk-2](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-unnamed-chunk-23.png) 




```r
p1 + geom_line(aes(time, escapement, group = reps), alpha = 0.1, col="darkgrey")
```

![plot of chunk unnamed-chunk-3](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-unnamed-chunk-33.png) 


### Computing additional statistics about the data


```r
crashed <- dt[time==OptTime, fishstock == 0, by=reps]
rewarded <- dt[time==OptTime, fishstock > xT, by=reps]
```






```r
dt <- data.table(dt, id=1:dim(dt)[1])
profits <- dt[, profit(fishstock, harvest), by=id]
```






```r
setkey(dt, id)
setkey(profits, id)
dt <- dt[profits]
setnames(dt, "V1", "profits")
setkey(dt, reps)
```






```r
total_profit <- dt[,sum(profits), by=reps]
total_profit <- total_profit + rewarded$V1 * reward 
```






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


```r
stats <- dt[ , mean_sdl(profits), by = time]
```



```
Error: columns of j don't evaluate to consistent types for each group: result for group 51 has column 1 type 'logical' but expecting type 'numeric'
```



```r
p1 + geom_line(dat=stats, aes(x=time, y=y), col="lightgrey") + 
  geom_ribbon(aes(x = time, ymin = ymin, ymax = ymax),
              fill = "darkred", alpha = 0.2, dat=stats)
```



```
Error: object 'stats' not found
```






```r
ggplot(dt, aes(total.profit, fill=crashed)) + geom_histogram(alpha=.8)
```

![plot of chunk unnamed-chunk-10](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-unnamed-chunk-103.png) 



