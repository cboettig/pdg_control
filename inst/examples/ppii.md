




# perfect policy, imperfect implementation 
Here's a trivial example I had in mind.  Trivial in the sense that the optimal policy is calculated assuming there will not be implementation error, but the reality has implementation error.  Because the model has alternate stable states, the errors drive the population over the threshold and cause it to crash.

I compare it against the same policy in which I tweak to be more conservative, and hence non-optimal.  Like the previous case, this example doesn't know about implementation error either, but it results in fewer crashes due to the non-optimal policy.

The optimal policy calculations are just straight-froward implementations of Reed, W.J., 1979, in which population growth is stochastic.

[*Reed 1979.  Optimal Escapement Levels in Stochastic and Deterministic Harvesting Models. Journal of Environmental Economics and Management. 6: 350-363.*]()



### Model setup 
Clear the workspace and load package dependencies. `pdgControl` is my implementation of these optimization methods and should be [installed from this repository](https://github.com/cboettig/pdg_control), see README. 



Define basic parameters of the economic optimal control problem.   We have a discrete economic discounting rate, and will solve the dynamic problem over a time window of 50 years.  In the discrete implementation we do not inforce the boundary condition, but rather put a value on meeting it.  The optimal solution may choose to not statisfy this boundary condition if the benefit outways this lost reward. 


```r
delta <- 0.1      # economic discounting rate
OptTime <- 50     # stopping time
reward <- 1       # bonus for satisfying the boundary condition
```




Use log-normal noise functions for growth noise, measurement error in the stock assessment, and implementation error, following the notation and definitions in [Sethi et al. (2005)](http://dx.doi.org/10.1016/j.jeem.2004.11.005).  To begin, we will allow only noise in growth, as in Reed 1979. 



```r
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
```




Chose the state equation / population dynamics function to have alternate stable states.  This is a Beverton-Holt like model with an Allee effect, a model based off of [Myers et al. (1995)](http://dx.doi.org/10.1126/science.269.5227.1106).  Note here we're just loading the function from the package.  The equilibrium value K is calculated explicitly from the model given this choice of parameters, as is the allee threshold.  We'll use the allee threshold as the final value `xT`. We'll start the model at the unharvested stochastic equilbrium size. 


```r
f <- Myer_harvest
pars <- c(1, 2, 6) 
p <- pars # shorthand 
K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
xT <- p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 # allee threshold
x0 <- K - sigma_g ^ 2 / 2 
```





We define a profit function with no stock effect for simplicity.  Profit is just linear in price, with some operating cost (which prevents strategies that put out more fishing effort than required when trying to catch all fish.)


```r
profit <- profit_harvest(price_fish = 1, 
                         cost_stock_effect = 0,
                         operating_cost = 0.1)
```




We solve the system numerically on a discrete grid. We'll consider all possible fish stocks between zero and twice the carrying capacity, and we'll consider harvest levels at the same resolution. 


```r
gridsize <- 100   # gridsize (discretized population)
x_grid <- seq(0, 2 * K, length = gridsize)  
h_grid <- x_grid  
```





### The perfect policy 
Having defined the problem, we are now ready to calculate the optimal policy by Bellman's stochastic dynamic programming solution.  We compute the stochastic transition matrices giving the probability that the stock goes from any value on the grid `x` at time `t` to any other value `y` at time `t+1`, for each possible harvest value.  Then we use this to compute the optimal harvest rate at each time interval in the system, solving backwards by dynamic programming.  The functions to do this are implemented in this package.


```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
```




### The imperfect implementation
Here we see how this policy performs over 100 replicates when implemented imperfectly.  The noise in the implementation was not part of the optimization.


```r
sigma_i <- 0.4 
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0, opt$D, z_g, z_m, z_i)
})
```




#### Outcome 
We summarize the results of the simulation in a tidy data table, ready for plotting.


```r
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps") # names are nice
```




This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown (solid line), along with the dotted line showing the allee threshold, below which the stock will go to zero (unless rescued stochastically). Colored dots behind the dynamics show the Optimal policy solution itself. 


```r
policy <- melt(opt$D)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
p6 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p6 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2, data=dt)
```

![plot of chunk fishstock_policy](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-fishstock_policy7.png) 


Calculate which of the replicates have crashed, which gives us a crude estimate of how badly the population has done. 


```r
crashed <- dt[time==as.integer(OptTime-1), fishstock < xT/4, by=reps]
```



A total of `40` crash.


Here we look at the dynamics of a single replicate, comparing the harvest and stock conditions.  Implementation errors give rise to violations of Reed's optimum escapement policy. 


```r
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, fishstock)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_line(aes(time, harvest), col="darkgreen") 
```

![plot of chunk rep](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-rep1.png) 


We can judge the performance on its own terms by looking at the distribution of total profit accrued in each simulation. Currently the simulation routine doesn't compute this on the fly, so it takes a bit of bookkeeping. 


```r
rewarded <- dt[time==OptTime, fishstock > xT, by=reps]

dt <- data.table(dt, id=1:dim(dt)[1])
profits <- dt[, profit(fishstock, harvest), by=id]

setkey(dt, id)
setkey(profits, id)
dt <- dt[profits]
setnames(dt, "V1", "profits")
setkey(dt, reps)

total_profit <- dt[,sum(profits), by=reps]
total_profit <- total_profit + rewarded$V1 * reward 

setkey(total_profit, reps)
setkey(crashed, reps)
setkey(rewarded, reps)
dt <- dt[total_profit]
dt <- dt[crashed]
dt <- dt[rewarded]
setnames(dt, c("V1", "V1.1", "V1.2"), c("total.profit", "crashed", "rewarded"))

##  totals
ggplot(dt, aes(total.profit, fill=crashed)) + geom_histogram(alpha=.8)
```

![plot of chunk profit](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-profit.png) 



## A non-optimal policy 
Let's adjust the optimal policy by a rule-of-thumb buffer, resulting in a non-optimal policy.


```r
buffer <- 0.05
safe_policy <- matrix(sapply(opt$D - buffer * length(h_grid), function(x) max(1, x)), ncol=dim(opt$D)[2])
```




This adds a `5` % buffer below the optimal harvest rate. 




```r
sims <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0, safe_policy, z_g, z_m, z_i)
})
```






```r
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "reps")
```






```r
policy <- melt(safe_policy)
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
p5 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)
p5 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2, data=dt)
```

![plot of chunk fishstock_policy2](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-fishstock_policy24.png) 




```r
crashed <- dt[time==as.integer(OptTime-1), fishstock < xT/4, by=reps]
```



A total of `8` crash.

Single replicate


```r
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, fishstock)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_line(aes(time, harvest), col="darkgreen") 
```

![plot of chunk rep2](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-rep21.png) 


We can look at the distribution of profits from this policy as well,


```r
rewarded <- dt[time==OptTime, fishstock > xT, by=reps]

dt <- data.table(dt, id=1:dim(dt)[1])
profits <- dt[, profit(fishstock, harvest), by=id]

setkey(dt, id)
setkey(profits, id)
dt <- dt[profits]
setnames(dt, "V1", "profits")
setkey(dt, reps)

total_profit <- dt[,sum(profits), by=reps]
total_profit <- total_profit + rewarded$V1 * reward 

setkey(total_profit, reps)
setkey(crashed, reps)
setkey(rewarded, reps)
dt <- dt[total_profit]
dt <- dt[crashed]
dt <- dt[rewarded]
setnames(dt, c("V1", "V1.1", "V1.2"), c("total.profit", "crashed", "rewarded"))

##  totals
ggplot(dt, aes(total.profit, fill=crashed)) + geom_histogram(alpha=.8)
```

![plot of chunk unnamed-chunk-2](http://www.carlboettiger.info/wp-content/uploads/2012/03/wpid-unnamed-chunk-213.png) 
