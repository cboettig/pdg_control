







```r
profit <- profit_harvest(price = 10, c0 = 30, c1 = 0)
```





```r
seed <- 123                 # Random seed (replicable results)
delta <- 0.01               # economic discounting rate
OptTime <- 100              # stopping time
gridsize <- 100              # grid size for fish stock and harvest rate (discretized population)
sigma_g <- 0.2              # Noise in population growth
reward <- 0                 # bonus for satisfying the boundary condition
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() 1         # No measurement noise, 
z_i <- function() 1         # No implemenation noise
f <- BevHolt                # Select the state equation
pars <- c(1.5, 0.05)        # parameters for the state equation
K <- (pars[1] - 1)/pars[2]  # Carrying capacity (for reference 
xT <- 0                     # boundary conditions
x0 <- K
x_grid <- seq(0.01, 1.2 * K, length = gridsize)  
h_grid <- seq(0.01, 0.8 * K, length = gridsize)  
```




```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
```





```r
L1 <- function(c2) function(h, h_prev)  c2 * abs(h - h_prev) 
free_increase <- function(c2) function(h, h_prev)  c2 * abs(min(h - h_prev, 0)) # increasing harvest is free
free_decrease <- function(c2) function(h, h_prev)  c2 * max(h - h_prev, 0) # decreasing harvest is free
fixed <-  function(c2) function(h, h_prev) c2 * as.numeric( !(h == h_prev) )
L2 <- function(c2) function(h, h_prev)  c2 * (h - h_prev) ^ 2
none <- function(h, h_prev)  0
penaltyfns <- list(L2=L2, L1=L1, free_decrease=free_decrease, fixed=fixed, free_increase=free_increase)
```




```r
apples <- c(L2=1.034, L1=1.724, free_decrease=1.724, fixed=8.276, free_increase=1.724)
xtable(apples)
```

```
Error: could not find function "xtable"
```




## Results

Solve the policy cost for the specified penalty function


```r
L2_policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                    profit, delta, reward, penalty = L2(apples["L2"]))
```



```r
fixed_policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
    profit, delta, reward, penalty = fixed(apples["fixed"]))
```



```r
L1_policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, reward, penalty = L1(apples["L1"]))
```



```r
free_increase_policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, reward, penalty =  free_increase(apples["free_increase"]))
```



```r
free_decrease_policy <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, profit, delta, reward, penalty = free_decrease(apples["free_decrease"]))
```


We also compare to the case in which costs of harvesting increase quadratically with effort; a common approach to create smoother policies.  


```r
quad_profit <- profit_harvest(price = 10, c0 = 30, c1 = apples["quad"]) 
quad_costs <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, quad_profit, delta, reward, penalty =  none)
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Warning: no non-missing arguments to max; returning -Inf
```

```
Error: incorrect number of dimensions
```




```r
sims <- list(
  L1 = simulate_optim(f, pars, x_grid, h_grid, x0, 
                      L1_policy$D, z_g, z_m, z_i, 
                      opt$D, profit=profit, penalty=L1(apples["L1"]), seed=seed), 
  L2 = simulate_optim(f, pars, x_grid, h_grid, x0, 
                      L2_policy$D, z_g, z_m, z_i, 
                      opt$D, profit=profit, penalty=L2(apples["L2"]), seed=seed),
  fixed = simulate_optim(f, pars, x_grid, h_grid, x0, 
                         fixed_policy$D, z_g, z_m, z_i, 
                         opt$D, profit=profit, penalty=fixed(apples["fixed"]), seed=seed),
  increase = simulate_optim(f, pars, x_grid, h_grid, x0, 
                            free_increase_policy$D, z_g, z_m, z_i, 
                            opt$D, profit=profit, penalty= free_increase(apples["increase"]), seed=seed),
  decrease = simulate_optim(f, pars, x_grid, h_grid, x0, 
                            free_decrease_policy$D, z_g, z_m, z_i, 
                            opt$D, profit=profit, penalty= free_decrease(apples["decrease"]), seed=seed),
  quad = simulate_optim(f, pars, x_grid, h_grid, x0, 
                            free_decrease_policy$D, z_g, z_m, z_i, 
                            opt$D, profit=quad_profit, penalty= none, seed=seed)
)
```




```r
#Make data tidy (melt), fast (data.tables), and nicely labeled.
dat <- melt(sims, id=names(sims[[1]]))  
dt <- data.table(dat)
setnames(dt, "L1", "penalty_fn") # names are nice
```


## Plots 




```r
p0 <- ggplot(dt) +
  geom_line(aes(time, harvest_alt), col="grey20", lwd=1) +
  geom_line(aes(time, harvest, col=penalty_fn, lty=penalty_fn))+ 
  labs(x="time", y="stock size", title = "Stock Dynamics")
p0
```

![plot of chunk p0](http://carlboettiger.info/assets/figures/2012-12-03-f4cdec09eb-p0.png) 



```r
p1 <- ggplot(dt) +
  geom_line(aes(time, alternate), col="grey20", lwd=1) +
  geom_line(aes(time, fishstock), col=rgb(0,0,1,.8)) + facet_wrap(~penalty_fn) + 
  labs(x="time", y="stock size", title = "Stock Dynamics")
p1
```

![plot of chunk p1](http://carlboettiger.info/assets/figures/2012-12-03-f4cdec09eb-p1.png) 



```r
p2 <- ggplot(dt) +
  geom_line(aes(time, harvest_alt), col="grey20", lwd=1)  +
  geom_line(aes(time, harvest), col=rgb(0,0,1,.8)) + 
  facet_wrap(~penalty_fn) + 
  labs(x="time", y="havest intensity (fish taken)", title = "Harvest Policy Dynamics")
p2
```

![plot of chunk p2](http://carlboettiger.info/assets/figures/2012-12-03-f4cdec09eb-p2.png) 


