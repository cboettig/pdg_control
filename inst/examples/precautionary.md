




# Precautionary Principle


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
x0 <- K - sigma_g ^ 2 / 2 
profit <- profit_harvest(price_fish = 1, cost_stock_effect = 0,
 operating_cost = 0.1)
x_grid <- seq(0, 2 * K, length = gridsize)  
h_grid <- x_grid  
```





### Calculate the transition matrix (with noise in growth only)      

Here we assume a precautionary high value of the threshold parameter.


```r
SDP_Mat <- determine_SDP_matrix(f, c(1,4,2), x_grid, h_grid, sigma_g )
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
```




### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above, with a threshold value that is actually lower than we assumed. 


```r
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, c(1,K,1), x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})
```




## Summarize and plot the results                                                   


```r
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
```




### Plots 
This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again, along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). 


```r
policy <- melt(opt$D)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
p6 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p6 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2, data=dt)
```

![plot of chunk fishstock_policy](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-fishstock_policy.png) 


Calculate which crashed


```r
crashed <- dt[time==as.integer(OptTime-1), fishstock < xT/4, by=reps]
```



A total of `5` crash.



### So is this more robust?
Here's the policy under higher intrinsic noise than it was designed for:


```r
sigma_g <- .25
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, c(1,K,1), x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})
```




Make data tidy (melt), fast (data.tables), and nicely labeled.


```r
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
```






```r
policy <- melt(opt$D)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
p6 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p6 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2, data=dt)
```

![plot of chunk unnamed-chunk-2](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-unnamed-chunk-211.png) 




```r
crashed <- dt[time==as.integer(OptTime-1), fishstock < xT/4, by=reps]
```



A total of `14` crash.


## Compare to the optimal policy performance when faced with additional growth noise


```r
pars <- c(1,K,1)
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g=.2)
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
sigma_g <- .25
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
crashed <- dt[time==as.integer(OptTime-1), fishstock == 0, by=reps]
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) + 
  geom_abline(intercept=xT, slope = 0, lty=2) 
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2)
```

![plot of chunk unnamed-chunk-4](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-unnamed-chunk-43.png) 


A total of `NaN &times; 10<sup>-Inf</sup>` crash.

