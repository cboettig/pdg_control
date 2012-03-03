





# Policy Costs 
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0





This example illustrates the impact of adding a cost to changing the harvest level between years 

### Define all parameters 


```r
delta <- 0.1      # economic discounting rate
OptTime <- 50     # stopping time
gridsize <- 100   # gridsize (discretized population)
sigma_g <- 0.2    # Noise in population growth
sigma_m <- 0.     # noise in stock assessment measurement
sigma_i <- 0.     # noise in implementation of the quota
reward <- 1       # bonus for satisfying the boundary condition
```




we'll use log normal noise functions


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
e_star <- (p[1] * sqrt(p[3]) - 2) / 2 ## Bifurcation point 
control <- "harvest"          # control variable can be harvest or effort 
price <- 1
cost <- .01
```




Our initial condition is the equilibrium size (note the stochastic deflation of mean)


```r
x0 <- K - sigma_g ^ 2 / 2 
```




and we use a harvest-based profit function with default parameters


```r
profit <-
function(x_grid, h_i){
  price_fish <- price
  stock_effect <- 0
  operating_cost <- cost * price_fish # trying to harvest more than costs
  sapply(x_grid, function(x_i){
    price_fish * min(h_i, x_i) - operating_cost * h_i  - stock_effect / x_i
  })
}
```




Set up the discrete grids for stock size and havest levels (which will use same resolution as for stock). 


```r
x_grid <- seq(0, 2 * K, length = gridsize)  
h_grid <- x_grid  
```




### Calculate the stochastic transition matrix
We calculate the stochastic transition matrix for the probability of going from any state \(x_t \) to any other state \(x_{t+1}\) the following year, for each possible choice of harvest \( h_t \).  This provides a look-up table for the dynamic programming calculations. Note that this only includes uncertainty in the growth rate (projected stock next year). 


```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
```




```r
require(snowfall) 
sfInit(parallel=TRUE, cpu=4)
SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps=999)
```



### Find the optimum by dynamic programming 
We use Bellman's dynamic programming algorithm to compute the optimal solution for all possible trajectories, ignoring potential policy costs as before.  We will later use this solution to compare against the optimal solution with policy costs.


```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                     profit, delta, reward=reward)
```



```
Error: incorrect number of dimensions
```




A modified algorithm lets us include a penalty of magnitude `P` and a functional form that can be an `L1` norm, `L2`  norm, `asymmetric` L1 norm, fixed cost, or `none` (no cost).  Here is an asymmetric norm example.  Note that this calculation is considerably slower. 


```r
policycost <- optim_policy(SDP_Mat, x_grid, h_grid, OptTime, xT, 
                    profit, delta, reward=reward, P=.3, penalty="asym")
```



```
Error: incorrect number of dimensions
```





### Simulate 
Now we'll simulate 100 replicates of this stochastic process under the optimal harvest policy determined above.  We use a modified simulation function that can simulate an alternate policy (the Reed optimum, where policy costs are zero, `opt$D` ) and a focal policy, `policycost$D`



```r
sims <- lapply(1:100, function(i)
  simulate_optim(f, pars, x_grid, h_grid, x0, policycost$D, z_g, z_m, z_i, opt$D)
  )
```



```
Error: object 'policycost' not found
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

Compare the optimal policy that involves this cost:


```r
policy <- melt(policycost$D)
```



```
Error: object 'policycost' not found
```



```r
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$fishstock) )
```



```
Error: object 'policy' not found
```



```r
p5 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=xT, slope=0, lty=2)
```



```
Error: object 'policy_zoom' not found
```



```r
p5 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1, data=dt)
```



```
Error: object 'p5' not found
```




Against the policy with no cost: 


```r
policy <- melt(opt$D)
```



```
Error: object 'opt' not found
```



```r
policy_zoom <- subset(policy, x_grid[Var1] < max(dt$alternate) )
```



```
Error: object 'policy' not found
```



```r
p6 <- ggplot(policy_zoom) + 
  geom_point(aes(Var2, (x_grid[Var1]), col=x_grid[Var1] - h_grid[value])) + 
  labs(x = "time", y = "fishstock") +
  scale_colour_gradientn(colours = rainbow(4)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_abline(intercept=xT, slope=0, lty=2)  
```



```
Error: object 'policy_zoom' not found
```



```r
p6 + geom_line(aes(time, alternate, group = reps), alpha = 0.1, data=dt)
```



```
Error: object 'p6' not found
```





Compare dynamics on a single replicate to see how this policy differs from the Reed policy. 


```r
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, fishstock)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_line(aes(time, harvest), col="darkgreen") 
```



```
Error: object 'reps' not found
```






