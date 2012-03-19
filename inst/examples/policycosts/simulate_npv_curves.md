




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






```r
delta <- 0.05     # economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
reward <- 0       # bonus for satisfying the boundary condition
```








```r
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
```








```r
f <- BevHolt                # Select the state equation
pars <- c(1.5, 0.05)             # parameters for the state equation
K <- (pars[1] - 1)/pars[2]  # Carrying capacity (for reference 
xT <- 0                     # boundary conditions
x0 <- K
```






```r
profit <- profit_harvest(price = 10, c0 = 30, c1 = 10)
```






```r
x_grid <- seq(0.01, 1.2 * K, length = gridsize)  
h_grid <- seq(0.01, 0.8 * K, length = gridsize)  
```






```r
L1 <- function(c2) function(h, h_prev)  c2 * abs(h - h_prev) 
asymmetric <- function(c2) function(h, h_prev)  c2 * max(h - h_prev, 0)
fixed <-  function(c2) function(h, h_prev) c2
L2 <- function(c2) function(h, h_prev)  c2 * (h - h_prev) ^ 2
```





## Loop over each penalty function with increasing penalty



```r
penaltyfns <- list(L2=L2, L1=L1, asy=asymmetric, fixed=fixed)
policies <- 
lapply(penaltyfns, function(penalty){
  c2 <- seq(0.01, 3, length.out = 15)
  policies <- 
  lapply(c2, function(c2){
    SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
    opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                         profit, delta, reward=reward)
    policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                        profit, delta, reward, penalty = penalty(c2))
    sims <- lapply(1:100, function(i)
      simulate_optim(f, pars, x_grid, h_grid, x0, policycost$D, z_g, z_m, z_i, opt$D, profit=profit, penalty=penalty(c2))
      )
    dat <- melt(sims, id=names(sims[[1]])) 
    names(dat) <- c(names(sims[[1]]), "reps")
    dat
  })

  dat <- melt(policies, id = names(policies[[1]]))
  dat$L1 <- c2[dat$L1]
  names(dat) <- c(names(policies[[1]]), "c2")
 dat
})
dat <- melt(policies, id = names(policies[[1]]))
names(dat) <- c(names(policies[[1]]), "penalty")
c2 <-  seq(0.01, 3, length.out = 15)
c2_index <- sapply(dat$c2, function(c) which(c2 %in% c))
dat2 <- data.frame(dat, c2_index = c2_index)
dt <- data.table(dat2)
setkey(dt,penalty,c2_index,reps)
```




## Plots

Variance trend


```r
vartrend <- 
sapply(c("L2", "L1", "fixed", "asy"), function(plty){
  sapply(1:15, function(i){
    pl <- quote(plty)
    mean(dt[penalty==eval(pl) & c2_index ==i,
         sd(harvest)/sd(fishstock), by=reps]$V1)
  })
})
df <- melt(data.frame(c2=c2, vartrend),id="c2") 
p1 <- ggplot(df) + geom_point(aes(c2, value, color=variable))+ geom_line(aes(c2, value, color=variable))
p1
```

![plot of chunk unnamed-chunk-2](http://farm8.staticflickr.com/7224/6850042286_ef81b74acc_o.png) 


Correlation trend


```r
cortrend <- sapply(c("L2", "L1", "fixed", "asy"), function(plty){
  sapply(1:15, function(i){
    pl <- quote(plty)
    mean(dt[penalty==eval(pl) & c2_index ==i,
         cor(harvest, fishstock), by=reps]$V1)
  })
})
df <- melt(data.frame(c2=c2, cortrend),id="c2") 
p2 <- ggplot(df) + geom_point(aes(c2, value, color=variable)) + geom_line(aes(c2, value, color=variable))
p2
```

![plot of chunk unnamed-chunk-3](http://farm8.staticflickr.com/7248/6996165783_41c9894bdb_o.png) 



Net present value plots


```r
discount <- (1-delta)^seq(1:50)
npv <- sapply(c("L2", "L1", "fixed", "asy"), function(plty){
  sapply(1:15, function(i){
    pl <- quote(plty)
    mean(dt[penalty==eval(pl) & c2_index ==i,
         sum(profit_fishing*discount), by=reps]$V1)
  })
})
df <- melt(data.frame(c2=c2, npv),id="c2") 
NPV0 <- mean(npv[1,])
p2 <- ggplot(df) + geom_point(aes(c2, (NPV0-value)/NPV0, color=variable)) + geom_line(aes(c2, (NPV0-value)/NPV0, color=variable))
p2
```

![plot of chunk unnamed-chunk-4](http://farm7.staticflickr.com/6042/6996166285_19e365b92a_o.png) 








s


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
delta <- 0.05     # economic discounting rate
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
pars <- c(1.5, 0.05)             # parameters for the state equation
K <- (pars[1] - 1)/pars[2]  # Carrying capacity (for reference 
xT <- 0                     # boundary conditions
x0 <- K
```




and we use a harvest-based profit function with default parameters


```r
profit <- profit_harvest(price = 10, c0 = 30, c1 = 10)
```




Set up the discrete grids for stock size and havest levels


```r
x_grid <- seq(0.01, 1.2 * K, length = gridsize)  
h_grid <- seq(0.01, 0.8 * K, length = gridsize)  
```





### Calculate the stochastic transition matrix
We calculate the stochastic transition matrix for the probability of going from any state \(x_t \) to any other state \(x_{t+1}\) the following year, for each possible choice of harvest \( h_t \).  This provides a look-up table for the dynamic programming calculations. Note that this only includes uncertainty in the growth rate (projected stock next year). 


```r
    SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
```



### Find the optimum by dynamic programming 

I've updated the algorithm to allow an arbitrary penalty function. Must be a function of the harvest and previous harvest. 


```r
L1 <- function(c2) function(h, h_prev)  c2 * abs(h - h_prev) 
policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                    profit, delta, reward, penalty = L1(.5))
```





### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above.  We use a modified simulation function that can simulate an alternate policy (the Reed optimum, where policy costs are zero, `opt$D` ) and a focal policy, `policycost$D`



```r
sims <- lapply(1:100, function(i)
  simulate_optim(f, pars, x_grid, h_grid, x0, policycost$D, z_g, z_m, z_i, opt$D)
  )
```



```
Error: argument "profit" is missing, with no default
```





Make data tidy (melt), fast (data.tables), and nicely labeled.



### Plots 

A single replicate, alternate dynamics should show the Reed optimum, while harvest/fishstock should show the impact of having policy costs. 


```r
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, alternate)) +
  geom_line(aes(time, fishstock), col="darkblue") +
  geom_line(aes(time, harvest), col="purple") + 
  geom_line(aes(time, harvest_alt), col="darkgreen") 
```



```
Error: object 'reps' not found
```





We can visualize the equilibrium policy for each possible harvest:



```r
policy <- sapply(1:length(h_grid), function(i) policycost$D[[i]][,1])
ggplot(melt(policy)) + 
  geom_point(aes(h_grid[Var2], (x_grid[Var1]), col=h_grid[value]-h_grid[Var2])) + 
    labs(x = "prev harvest", y = "fishstock") +
      scale_colour_gradientn(colours = rainbow(4)) 
```

![plot of chunk unnamed-chunk-5](http://farm7.staticflickr.com/6239/6850051300_8c3311f17c_o.png) 


Here we plot previous harvest against the recommended harvest, coloring by stocksize.  Note this swaps the y axis from above with the color density.  Hence each x-axis value has all possible colors, but they map down onto a subset of optimal harvest values (depending on their stock). 


```r
policy <- sapply(1:length(h_grid), function(i) policycost$D[[i]][,1])
ggplot(melt(policy)) + 
  geom_point(aes(h_grid[Var2], (h_grid[value]), col = x_grid[Var1]), position=position_jitter(w=.005,h=.005), alpha=.5) + 
    labs(x = "prev harvest", y = "harvest") +
      scale_colour_gradientn(colours = rainbow(4)) 
```

![plot of chunk unnamed-chunk-6](http://farm8.staticflickr.com/7239/6996174525_e71b4ae557_o.png) 



### Profits


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
Error: object 'fishstock' not found
```




Merge in profits to data.table (should be a way to avoid having to do these joins?)


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




merge in total profits to data.table


```r
total_profit <- dt[,sum(profits), by=reps]
```



```
Error: object 'profits' not found
```



```r
setkey(total_profit, reps)
```



```
Error: object 'total_profit' not found
```



```r
setkey(dt, reps)
```



```
Error: x is not a data.table
```



```r
dt <- dt[total_profit]
```



```
Error: object 'total_profit' not found
```



```r
setnames(dt, "V1", "total.profit")
```



```
Error: x is not a data.table
```






```r
ggplot(dt, aes(total.profit)) + geom_histogram(alpha=.8)
```



```
Error: ggplot2 doesn't know how to deal with data of class function
```






```r
save(list=ls(), file="L1.rda")
```




The mean dynamics of the state


```r
stats <- dt[ , mean_sdl(fishstock), by = time]
```



```
Error: object 'fishstock' not found
```



```r
ggplot(stats) +   geom_ribbon(aes(x = time, ymin = ymin, ymax = ymax),
                fill = "darkblue", alpha = 0.2, dat=stats) +
                geom_line(aes(x=time, y=y), lwd=1) 
```



```
Error: object 'stats' not found
```




The mean dynamics of the control


```r
stats <- dt[ , mean_sdl(harvest), by = time]
```



```
Error: object 'harvest' not found
```



```r
ggplot(stats) +  geom_ribbon(aes(x = time, ymin = ymin, ymax = ymax),
                fill = "darkblue", alpha = 0.2) +
                geom_line(aes(x=time, y=y), lwd=1) 
```



```
Error: object 'stats' not found
```


