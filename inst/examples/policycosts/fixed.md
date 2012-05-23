





# Policy Costs fixed fee 
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

## Setup the system



```r
rm(list = ls())
require(pdgControl)
require(reshape2)
require(ggplot2)
require(data.table)

delta <- 0.05  # economic discounting rate
OptTime <- 50  # stopping time
gridsize <- 50  # gridsize (discretized population)
sigma_g <- 0.2  # Noise in population growth
reward <- 0  # bonus for satisfying the boundary condition

z_g <- function() rlnorm(1, 0, sigma_g)  # mean 1
z_m <- function() 1
z_i <- function() 1

f <- BevHolt  # Select the state equation
pars <- c(1.5, 0.05)  # parameters for the state equation
K <- (pars[1] - 1)/pars[2]  # Carrying capacity (for reference
xT <- 0  # boundary conditions
x0 <- K


profit <- profit_harvest(price = 10, c0 = 30, c1 = 0)

x_grid <- seq(0.01, 1.2 * K, length = gridsize)
h_grid <- seq(0.01, 0.8 * K, length = gridsize)


SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g)
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, 
    delta, reward = reward)
```








```r
L1 <- function(c2) function(h, h_prev) c2 * abs(h - h_prev)
free_increase <- function(c2) function(h, h_prev) c2 * abs(min(h - 
    h_prev, 0))  # increasing harvest is free
free_decrease <- function(c2) function(h, h_prev) c2 * max(h - h_prev, 
    0)  # decreasing harvest is free
fixed <- function(c2) function(h, h_prev) c2 * as.numeric(!(h == 
    h_prev))
L2 <- function(c2) function(h, h_prev) c2 * (h - h_prev)^2
```




Solve the policy cost for the specified penalty function



```r
c2 <- 9.103
penalty <- fixed(c2)
policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
    profit, delta, reward, penalty = penalty)
cache = FALSE
```





### Simulate 

Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above.  We use a modified simulation function that can simulate an alternate policy (the Reed optimum, where policy costs are zero, `opt$D` ) and a focal policy, `policycost$D`



```r
sims <- lapply(1:100, function(i) simulate_optim(f, pars, x_grid, 
    h_grid, x0, policycost$D, z_g, z_m, z_i, opt$D, profit = profit, penalty = penalty))
```




Make data tidy (melt), fast (data.tables), and nicely labeled.



```r
dat <- melt(sims, id = names(sims[[1]]))
dt <- data.table(dat)
setnames(dt, "L1", "reps")  # names are nice
```




# Plots 

A single replicate, alternate dynamics show the Reed optimum, while harvest/fishstock show the impact of having policy costs. 



```r
ggplot(subset(dt, reps == 1)) + geom_line(aes(time, alternate)) + 
    geom_line(aes(time, fishstock), col = "darkblue") + geom_line(aes(time, 
    harvest), col = "purple") + geom_line(aes(time, harvest_alt), col = "darkgreen")
```

![plot of chunk rep1](http://farm8.staticflickr.com/7100/7253144964_6f2dcd6c93_o.png) 


A second replicate



```r
ggplot(subset(dt, reps == 2)) + geom_line(aes(time, alternate)) + 
    geom_line(aes(time, fishstock), col = "darkblue") + geom_line(aes(time, 
    harvest), col = "purple") + geom_line(aes(time, harvest_alt), col = "darkgreen")
```

![plot of chunk rep2](http://farm8.staticflickr.com/7212/7253145562_e9f9c2dbaa_o.png) 


## Profits 



```r
ggplot(subset(dt, reps == 1)) + geom_line(aes(time, profit_fishing)) + 
    geom_line(aes(time, policy_cost), col = "darkblue")
```

![plot of chunk rep1profit](http://farm9.staticflickr.com/8153/7253145890_d400249b8e_o.png) 




```r
costs <- dt[, sum(policy_cost), by = reps]
profits <- dt[, sum(profit_fishing), by = reps]

qplot(costs$V1)
```

![plot of chunk policycost](http://farm9.staticflickr.com/8153/7253146290_66b66c8d5a_o.png) 

```r
qplot(profits$V1)
```

![plot of chunk policycost](http://farm8.staticflickr.com/7091/7253146626_8469654821_o.png) 

```r
qplot(profits$V1 - costs$V1)
```

![plot of chunk policycost](http://farm8.staticflickr.com/7245/7253146924_92f0314d90_o.png) 




We can visualize the equilibrium policy for each possible harvest:



```r
policy <- sapply(1:length(h_grid), function(i) policycost$D[[i]][, 
    1])
ggplot(melt(policy)) + geom_point(aes(h_grid[Var2], (x_grid[Var1]), 
    col = h_grid[value] - h_grid[Var2])) + labs(x = "prev harvest", y = "fishstock") + 
    scale_colour_gradientn(colours = rainbow(4))
```

![plot of chunk policy](http://farm8.staticflickr.com/7093/7253147292_c1975b2b51_o.png) 


Here we plot previous harvest against the recommended harvest, coloring by stocksize.  Note this swaps the y axis from above with the color density.  Hence each x-axis value has all possible colors, but they map down onto a subset of optimal harvest values (depending on their stock). 



```r
policy <- sapply(1:length(h_grid), function(i) policycost$D[[i]][, 
    1])
ggplot(melt(policy)) + geom_point(aes(h_grid[Var2], (h_grid[value]), 
    col = x_grid[Var1]), position = position_jitter(w = 0.005, h = 0.005), alpha = 0.5) + 
    labs(x = "prev harvest", y = "harvest") + scale_colour_gradientn(colours = rainbow(4))
```

![plot of chunk harvestchanges](http://farm9.staticflickr.com/8023/7253147656_53fa29f7b3_o.png) 


