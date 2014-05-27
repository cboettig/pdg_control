




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
f <- RickerAllee
K <- 4 
xT <- 1 # final value, also allee threshold
pars <- c(1, K, xT) 
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

![plot of chunk fishstock](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-fishstock23.png) 


### Computing additional statistics about the data
In this section we add some additional information to our data.table on the profits obtained by each replicate.  The algorithm has supposedly maximized the expected profit, so it is useful to look at both the mean total profit and the distribution.  Despite this maximization, the distribution can be rather lop-sided or even bimodal. 

Which replicates crashed?  Which met the boundary requirment and recieved the reward value at the end?


```r
crashed <- dt[time==as.integer(OptTime-1), fishstock < xT/4, by=reps]
rewarded <- dt[time==OptTime, fishstock > xT, by=reps]
```




A total of `21` crash.



## Compare to a non-optimal solution
Compare another model, that likewise assumes no implementation error, and also makes a mistake in its estimate of the growth parameter, making it conservative rather than optimal.




```r
sigma_i <- 0
SDP_Mat <- determine_SDP_matrix(f, c(1,4,2), x_grid, h_grid, sigma_g )
```






```r
nonopt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
```





### Simulate 
For the simulated implementation, we add the same implementation error back, and we restore biological allee threshold to it's true value. 


```r
sigma_i <- .4
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, c(1,4,1), x_grid, h_grid, x0, nonopt$D, z_g, z_m, z_i)
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

![plot of chunk unnamed-chunk-1](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-unnamed-chunk-114.png) 


### Computing additional statistics about the data


```r
crashed <- dt[time==as.integer(OptTime-1), fishstock < xT/4, by=reps]
rewarded <- dt[time==OptTime, fishstock > xT, by=reps]
```



A total of `22` crash.


