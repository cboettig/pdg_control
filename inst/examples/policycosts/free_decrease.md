






# Asymmetric Policy Costs 
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

## Setup the system



```r
source("setup.R")
```





```r
L1 <- function(c2) function(h, h_prev) c2 * abs(h - h_prev)
free_increase <- function(c2) function(h, h_prev) c2 * abs(min(h - 
    h_prev, 0))  # increasing harvest is free
free_decrease <- function(c2) function(h, h_prev) c2 * max(h - h_prev, 
    0)  # decreasing harvest is free
fixed <- function(c2) function(h, h_prev) c2
L2 <- function(c2) function(h, h_prev) c2 * (h - h_prev)^2
```




Solve the policy cost for the specified penalty function



```r
c2 <- 4
penalty <- free_decrease(c2)
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




### Plots 

A single replicate, alternate dynamics show the Reed optimum, while harvest/fishstock show the impact of having policy costs. 



```r
ggplot(subset(dt, reps == 1)) + geom_line(aes(time, alternate)) + 
    geom_line(aes(time, fishstock), col = "darkblue") + geom_line(aes(time, 
    harvest), col = "purple") + geom_line(aes(time, harvest_alt), col = "darkgreen")
```

![plot of chunk rep1](http://farm8.staticflickr.com/7197/7134149057_a3d16263e6_o.png) 


A second replicate



```r
ggplot(subset(dt, reps == 2)) + geom_line(aes(time, alternate)) + 
    geom_line(aes(time, fishstock), col = "darkblue") + geom_line(aes(time, 
    harvest), col = "purple") + geom_line(aes(time, harvest_alt), col = "darkgreen")
```

![plot of chunk rep2](http://farm8.staticflickr.com/7237/6988065216_f600d96f8d_o.png) 


We can visualize the equilibrium policy for each possible harvest:



```r
policy <- sapply(1:length(h_grid), function(i) policycost$D[[i]][, 
    1])
ggplot(melt(policy)) + geom_point(aes(h_grid[Var2], (x_grid[Var1]), 
    col = h_grid[value] - h_grid[Var2])) + labs(x = "prev harvest", y = "fishstock") + 
    scale_colour_gradientn(colours = rainbow(4))
```

![plot of chunk policy](http://farm9.staticflickr.com/8013/7134149517_8eccbb1248_o.png) 


Here we plot previous harvest against the recommended harvest, coloring by stocksize.  Note this swaps the y axis from above with the color density.  Hence each x-axis value has all possible colors, but they map down onto a subset of optimal harvest values (depending on their stock). 



```r
policy <- sapply(1:length(h_grid), function(i) policycost$D[[i]][, 
    1])
ggplot(melt(policy)) + geom_point(aes(h_grid[Var2], (h_grid[value]), 
    col = x_grid[Var1]), position = position_jitter(w = 0.005, h = 0.005), alpha = 0.5) + 
    labs(x = "prev harvest", y = "harvest") + scale_colour_gradientn(colours = rainbow(4))
```

![plot of chunk harvestchanges](http://farm8.staticflickr.com/7135/7134149855_ee35421052_o.png) 


## Profits



```r
save(list = ls(), file = "stochastic_norms.rda")
```



