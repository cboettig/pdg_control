




# Increasing costs  
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




### Define all parameters 


```r
delta <- 0.05     # economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
reward <- 0       # bonus for satisfying the boundary condition


c2 <- 1           # Policy adjustment cost
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
profit <- profit_harvest(price = 10, c0 = 30, c1 = 0.0)
```




Set up the discrete grids for stock size and havest levels


```r
x_grid <- seq(0, 1.2 * K, length = gridsize)  
h_grid <- seq(0, 0.8 * K, length = gridsize)  
```





Specify the functional form of the policy costs:


```r
L2 <- function(h, h_prev)  c2 * (h - h_prev) ^ 2
L1 <- function(h, h_prev)  c2 * abs(h - h_prev) 
asymmetric <- function(h, h_prev)  c2 * max(h_grid[h_prev] - h_grid[i], 0)
fixed <-  function(h, h_prev) c2
```




### Calculate the stochastic transition matrix
We calculate the stochastic transition matrix for the probability of going from any state \(x_t \) to any other state \(x_{t+1}\) the following year, for each possible choice of harvest \( h_t \).  This provides a look-up table for the dynamic programming calculations. Note that this only includes uncertainty in the growth rate (projected stock next year). 


```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
```



### Find the optimum by dynamic programming 
We use Bellman's dynamic programming algorithm to compute the optimal solution for all possible trajectories, ignoring potential policy costs as before.  We will later use this solution to compare against the optimal solution with policy costs.


```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
```




A modified algorithm lets us include a penalty of magnitude `P` and a functional form that can be an `L1` norm, `L2`  norm, `asymmetric` L1 norm (costly to lower harvest rates), fixed cost, or `none` (no cost).  Here is an asymmetric norm example.  Note that this calculation is considerably slower. 


```r
policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                    profit, delta, reward, penalty = L2)
```





### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above.  We use a modified simulation function that can simulate an alternate policy (the Reed optimum, where policy costs are zero, `opt$D` ) and a focal policy, `policycost$D`



```r
sims <- lapply(1:100, function(i)
  simulate_optim(f, pars, x_grid, h_grid, x0, policycost$D, z_g, z_m, z_i, opt$D, penalty=L2)
  )
```



```
Error: argument "profit" is missing, with no default
```





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

A single replicate, alternate dynamics should show the Reed optimum, while harvest/fishstock should show the impact of having policy costs. 


```r
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, alternate)) +
  geom_line(aes(time, fishstock), col="darkblue") +
  geom_abline(intercept=opt$S, slope = 0) +
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
  geom_point(aes(h_grid[Var2], (x_grid[Var1]), col=h_grid[value] - h_grid[Var2])) + 
    labs(x = "prev harvest", y = "fishstock") +
      scale_colour_gradientn(colours = rainbow(4)) 
```

![plot of chunk unnamed-chunk-2](http://farm8.staticflickr.com/7184/6837543660_35fce1a893_o.png) 


