


```r
rm(list = ls())
require(pdgControl)
require(reshape2)
require(ggplot2)
require(data.table)
```



```r
opts_knit$set(upload.fun = socialR::notebook.url)
opts_chunk$set(dev.args = list(bg = "transparent"))
theme_set(theme_bw())
theme_update(panel.background = element_rect(fill = "transparent", 
    colour = NA), plot.background = element_rect(fill = "transparent", colour = NA))
```




```r
delta <- 0.05  # economic discounting rate
OptTime <- 50  # stopping time
gridsize <- 50  # gridsize (discretized population)
sigma_g <- 0.2  # Noise in population growth
reward <- 0  # bonus for satisfying the boundary condition
z_g <- function() rlnorm(1, 0, sigma_g)  # mean 1
z_m <- function() 1  # measurement noise, none, but part of simulation routine
z_i <- function() 1
f <- BevHolt  # Select the state equation
pars <- c(1.5, 0.05)  # parameters for the state equation
K <- (pars[1] - 1)/pars[2]  # Carrying capacity (for reference
xT <- 0  # boundary conditions
x0 <- K
profit <- profit_harvest(price = 10, c0 = 0, c1 = 0)
x_grid <- seq(0.01, 1.2 * K, length = gridsize)
h_grid <- seq(0.01, 0.8 * K, length = gridsize)
```




```r
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
L2_policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, 
    delta, reward, penalty = L2(0.7586))
```




```r
fixed_policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
    profit, delta, reward, penalty = fixed(9.103))
```




```r
L1_policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, 
    delta, reward, penalty = L1(1.52))
```




```r
free_increase_policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, 
    xT, profit, delta, reward, penalty = free_increase(3.04))
```




```r
free_decrease_policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, 
    xT, profit, delta, reward, penalty = free_decrease(3.04))
```




```r
sims <- list(L1 = simulate_optim(f, pars, x_grid, h_grid, x0, L1_policy$D, 
    z_g, z_m, z_i, opt$D, profit = profit, penalty = L1(1.52), seed = 111), 
    L2 = simulate_optim(f, pars, x_grid, h_grid, x0, L2_policy$D, z_g, z_m, 
        z_i, opt$D, profit = profit, penalty = L2(0.7586), seed = 111), fixed = simulate_optim(f, 
        pars, x_grid, h_grid, x0, fixed_policy$D, z_g, z_m, z_i, opt$D, profit = profit, 
        penalty = fixed(9.103), seed = 111), increase = simulate_optim(f, pars, 
        x_grid, h_grid, x0, free_decrease_policy$D, z_g, z_m, z_i, opt$D, profit = profit, 
        penalty = free_increase(1.52), seed = 111), decrease = simulate_optim(f, 
        pars, x_grid, h_grid, x0, free_decrease_policy$D, z_g, z_m, z_i, opt$D, 
        profit = profit, penalty = free_decrease(1.52), seed = 111))
```


Make data tidy (melt), fast (data.tables), and nicely labeled.


```r
dat <- melt(sims, id = names(sims[[1]]))
dt <- data.table(dat)
setnames(dt, "L1", "penalty_fn")  # names are nice
```


# Plots 




```r
p0 <- ggplot(dt) + geom_line(aes(time, alternate), col = "grey20", 
    lwd = 1) + geom_line(aes(time, fishstock, col = penalty_fn, lty = penalty_fn)) + 
    labs(x = "time", y = "stock size", title = "Stock Dynamics")
p0
```

![plot of chunk unnamed-chunk-4](http://carlboettiger.info/assets/figures/2012-11-27-89f3cb4cfe-unnamed-chunk-4.png) 



```r
p1 <- ggplot(dt) + geom_line(aes(time, alternate), col = "grey20", 
    lwd = 1) + geom_line(aes(time, fishstock), col = rgb(0, 0, 1, 0.8)) + facet_wrap(~penalty_fn) + 
    labs(x = "time", y = "stock size", title = "Stock Dynamics")
p1
```

![plot of chunk unnamed-chunk-5](http://carlboettiger.info/assets/figures/2012-11-27-89f3cb4cfe-unnamed-chunk-5.png) 



```r
p2 <- ggplot(dt) + geom_line(aes(time, harvest_alt), col = "grey20", 
    lwd = 1) + geom_line(aes(time, harvest), col = rgb(0, 0, 1, 0.8)) + facet_wrap(~penalty_fn) + 
    labs(x = "time", y = "havest intensity (fish taken)", title = "Harvest Policy Dynamics")
p2
```

![plot of chunk unnamed-chunk-6](http://carlboettiger.info/assets/figures/2012-11-27-89f3cb4cfe-unnamed-chunk-6.png) 


