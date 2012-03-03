




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
 operating_cost = 0.1 * price)
```



```
Error: unused argument(s) (price_fish = 1, cost_stock_effect = 0, operating_cost = 0.1 * price)
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



```
Error: object 'profit' not found
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



```
Error: object 'opt' not found
```




## Summarize and plot the results                                                   
Make data tidy (melt), fast (data.tables), and nicely labeled.


```r
dat <- melt(sims, id=names(sims[[1]]))  
```



```
Error: object 'sims' not found
```



```r
dt <- data.table(dat)
```



```
Error: object 'dat' not found
```



```r
setnames(dt, "L1", "reps") # names are nice
```



```
Error: x is not a data.table
```




### Plots 
This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again, along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). 


```r
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) + 
  geom_abline(intercept=xT, slope = 0, lty=2) 
```



```
Error: ggplot2 doesn't know how to deal with data of class function
```



```r
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
```



```
Error: object 'p1' not found
```




We can also look at the harvest dynamics:


```r
p1 + geom_line(aes(time, harvest, group = reps), alpha = 0.1, col="darkgreen")
```



```
Error: object 'p1' not found
```




This strategy is supposed to be a constant-escapement strategy. We can visualize the escapement and see if it is less variable than fish stock, and if it is near Reed's S: 


```r
p1 + geom_line(aes(time, escapement, group = reps), alpha = 0.1, col="darkgrey")
```



```
Error: object 'p1' not found
```




### Computing additional statistics about the data
In this section we add some additional information to our data.table on the profits obtained by each replicate.  The algorithm has supposedly maximized the expected profit, so it is useful to look at both the mean total profit and the distribution.  Despite this maximization, the distribution can be rather lop-sided or even bimodal. 

Which replicates crashed?  Which met the boundary requirment and recieved the reward value at the end?


```r
crashed <- dt[time==OptTime, fishstock == 0, by=reps]
```



```
Error: comparison (1) is possible only for atomic and list types
```



```r
rewarded <- dt[time==OptTime, fishstock > xT, by=reps]
```



```
Error: comparison (1) is possible only for atomic and list types
```




A total of 


```r
sum(crashed$V1)
```



```
Error: object 'crashed' not found
```



crash.

Let's compute the profits at each time-step for each replicate. 
Using `data.table` to evaluate our profit function over the stock and harvest levels requires indexing our data:



```r
dt <- data.table(dt, id=1:dim(dt)[1])
```



```
Error: argument of length 0
```



```r
profits <- dt[, profit(fishstock, harvest), by=id]
```



```
Error: could not find function "profit"
```




Merging this calculation back into our data table using fast join (needs to define 'id' as a key on which to match things up though). 


```r
setkey(dt, id)
```



```
Error: x is not a data.table
```



```r
setkey(profits, id)
```



```
Error: object 'profits' not found
```



```r
dt <- dt[profits]
```



```
Error: object 'profits' not found
```



```r
setnames(dt, "V1", "profits")
```



```
Error: x is not a data.table
```



```r
setkey(dt, reps)
```



```
Error: x is not a data.table
```




Compute total profit by summing over each timeseries (including the reward for satisfying the terminal boundary condition, if any). 



```r
total_profit <- dt[,sum(profits), by=reps]
```



```
Error: object 'profits' not found
```



```r
total_profit <- total_profit + rewarded$V1 * reward 
```



```
Error: object 'total_profit' not found
```





Add these three columns to the data.table (fast join and re-label):


```r
setkey(total_profit, reps)
```



```
Error: object 'total_profit' not found
```



```r
setkey(crashed, reps)
```



```
Error: object 'crashed' not found
```



```r
setkey(rewarded, reps)
```



```
Error: object 'rewarded' not found
```



```r
dt <- dt[total_profit]
```



```
Error: object 'total_profit' not found
```



```r
dt <- dt[crashed]
```



```
Error: object 'crashed' not found
```



```r
dt <- dt[rewarded]
```



```
Error: object 'rewarded' not found
```



```r
setnames(dt, c("V1", "V1.1", "V1.2"), c("total.profit", "crashed", "rewarded"))
```



```
Error: x is not a data.table
```






#### Profit plots
Since the optimal strategy maximizes expected profit, it may be more useful to look at the distribution statistics of profit over time:


```r
stats <- dt[ , mean_sdl(profits), by = time]
```



```
Error: object 'profits' not found
```



```r
p1 + geom_line(dat=stats, aes(x=time, y=y), col="lightgrey") + 
  geom_ribbon(aes(x = time, ymin = ymin, ymax = ymax),
              fill = "darkred", alpha = 0.2, dat=stats)
```



```
Error: object 'p1' not found
```





Total profits


```r
ggplot(dt, aes(total.profit, fill=crashed)) + geom_histogram(alpha=.8)
```



```
Error: ggplot2 doesn't know how to deal with data of class function
```




## Compare to a non-optimal solution
Compare another model, that likewise assumes no implementation error, and also makes a mistake in its estimate of the growth parameter, making it conservative rather than optimal.



```r
pars <- c(1, K, 1.5)
```







```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
```



```
Error: object 'profit' not found
```




For the simulated implementation, we add the same implementation error back, and we restore biological growth noise to it's true value


```r
sigma_i <- 0.4
```





### Simulate 
Now we simulate as before


```r
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})
```



```
Error: object 'opt' not found
```




## Summarize and plot the results                                                  
Using the code above, recreate the plots for this policy and simulation: 


```r
dat <- melt(sims, id=names(sims[[1]]))  
```



```
Error: object 'sims' not found
```



```r
dt <- data.table(dat)
```



```
Error: object 'dat' not found
```



```r
setnames(dt, "L1", "reps") # names are nice
```



```
Error: x is not a data.table
```




### Plots 


```r
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) + 
  geom_abline(intercept=xT, slope = 0, lty=2) 
```



```
Error: ggplot2 doesn't know how to deal with data of class function
```



```r
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
```



```
Error: object 'p1' not found
```






```r
p1 + geom_line(aes(time, harvest, group = reps), alpha = 0.1, col="darkgreen")
```



```
Error: object 'p1' not found
```






```r
p1 + geom_line(aes(time, escapement, group = reps), alpha = 0.1, col="darkgrey")
```



```
Error: object 'p1' not found
```




### Computing additional statistics about the data


```r
crashed <- dt[time==OptTime, fishstock == 0, by=reps]
```



```
Error: comparison (1) is possible only for atomic and list types
```



```r
rewarded <- dt[time==OptTime, fishstock > xT, by=reps]
```



```
Error: comparison (1) is possible only for atomic and list types
```



A total of 


```r
sum(crashed$V1)
```



```
Error: object 'crashed' not found
```



crash.



```r
dt <- data.table(dt, id=1:dim(dt)[1])
```



```
Error: argument of length 0
```



```r
profits <- dt[, profit(fishstock, harvest), by=id]
```



```
Error: could not find function "profit"
```






```r
setkey(dt, id)
```



```
Error: x is not a data.table
```



```r
setkey(profits, id)
```



```
Error: object 'profits' not found
```



```r
dt <- dt[profits]
```



```
Error: object 'profits' not found
```



```r
setnames(dt, "V1", "profits")
```



```
Error: x is not a data.table
```



```r
setkey(dt, reps)
```



```
Error: x is not a data.table
```






```r
total_profit <- dt[,sum(profits), by=reps]
```



```
Error: object 'profits' not found
```



```r
total_profit <- total_profit + rewarded$V1 * reward 
```



```
Error: object 'total_profit' not found
```






```r
setkey(total_profit, reps)
```



```
Error: object 'total_profit' not found
```



```r
setkey(crashed, reps)
```



```
Error: object 'crashed' not found
```



```r
setkey(rewarded, reps)
```



```
Error: object 'rewarded' not found
```



```r
dt <- dt[total_profit]
```



```
Error: object 'total_profit' not found
```



```r
dt <- dt[crashed]
```



```
Error: object 'crashed' not found
```



```r
dt <- dt[rewarded]
```



```
Error: object 'rewarded' not found
```



```r
setnames(dt, c("V1", "V1.1", "V1.2"), c("total.profit", "crashed", "rewarded"))
```



```
Error: x is not a data.table
```




#### Profit plots


```r
stats <- dt[ , mean_sdl(profits), by = time]
```



```
Error: object 'profits' not found
```



```r
p1 + geom_line(dat=stats, aes(x=time, y=y), col="lightgrey") + 
  geom_ribbon(aes(x = time, ymin = ymin, ymax = ymax),
              fill = "darkred", alpha = 0.2, dat=stats)
```



```
Error: object 'p1' not found
```






```r
ggplot(dt, aes(total.profit, fill=crashed)) + geom_histogram(alpha=.8)
```



```
Error: ggplot2 doesn't know how to deal with data of class function
```





