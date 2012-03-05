




# perfect policy, imperfect implementation 
Compare to non-optimal, rule-of-thumb policy.

### Model setup 
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
e_star <- (p[1] * sqrt(p[3]) - 2) / 2 ## Bifurcation point, for reference 
x0 <- K - sigma_g ^ 2 / 2 
```






```r
profit <- profit_harvest(price_fish = 1, 
                         cost_stock_effect = 0,
                         operating_cost = 0.1)
```






```r
x_grid <- seq(0, 2 * K, length = gridsize)  
h_grid <- x_grid  
```





### The perfect policy 
Calculate the optimal policy


```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
```




### The imperfect implementation

Implementation the optimal polict with implementation noise 


```r
sigma_i <- 0.4 
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})
```




#### Outcome 


```r
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
```




This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown (solid line), along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). 


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

![plot of chunk fishstock_policy](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-fishstock_policy4.png) 


Calculate which crashed


```r
crashed <- dt[time==as.integer(OptTime-1), fishstock < xT/4, by=reps]
```



A total of `35` crash.



### A non-optimal policy 
Let's adjust the optimal policy by a rule-of-thumb buffer, resulting in a non-optimal policy.


```r
buffer <- 0.1
safe_policy <- matrix(sapply(opt$D - buffer * length(h_grid), function(x) max(0, x)), ncol=dim(opt$D)[2])
```




This adds a `10` % buffer below the optimal harvest rate. 




```r
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, c(1,K,1), x_grid, h_grid, x0, safe_policy, z_g, z_m, z_i)
})
```



```
Error: replacement has length zero
```






```r
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
```






```r
policy <- melt(safe_policy)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
p6 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p6 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2, data=dt)
```

![plot of chunk fishstock_policy2](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-fishstock_policy22.png) 




```r
crashed <- dt[time==as.integer(OptTime-1), fishstock < xT/4, by=reps]
```



A total of `35` crash.



