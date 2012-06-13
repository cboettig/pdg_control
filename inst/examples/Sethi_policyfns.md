




# Sethi Model
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

Implements a numerical version of the SDP described in (Sethi _et. al._ 2005).


Clear the workspace and load package dependencies: 




Chose the state equation / population dynamics function



```r
f <- BevHolt
```




That is, \\( f(x,h) = \frac{A x}{1 + B x} \\)

We must now define parameters for the function.  Note that the positive stationary root of the model is given by \\( \frac{A-1}{B} \\), which we'll store for future reference as `K`.  



```r
pars <- c(1.5, 0.05)
K <- (pars[1] - 1)/pars[2]
```




and we use a harvest-based profit function with default parameters



```r
profit <- profit_harvest(price=1, c0 = 0.01) 
```




The `profit_harvest` function has the form \\( \Pi = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x} \\), conditioned on \\( h > x \\) and \\(x > 0 \\). We set up the discrete grids for stock size and havest levels (which will use same resolution as for stock), in order to calculate the SDP solution.   Here we set the gridsize to 100.  



```r
x_grid <- seq(0, 1.5 * K, length = 100)  
h_grid <- x_grid  
```





## Scenarios

### Large Measurement error



```r
sigma_g <- 0.01    # Noise in population growth
sigma_m <- 0.25     # noise in stock assessment measurement
sigma_i <- 0.01     # noise in implementation of the quota
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
```






```r
require(snowfall) 
sfInit(parallel=TRUE, cpu=16)
```



```
R Version:  R version 2.14.1 (2011-12-22) 

```



```r
SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps=19999)
```



```
Library ggplot2 loaded.
```



Note that `SDP_Mat` is specified from the calculation above, as are our grids and our profit function. `OptTime` is the stopping time.  `xT` specifies a boundary condition at the stopping time. A reward for meeting this boundary must be specified for it to make any difference.  `delta` indicates the economic discount rate. Again, details are in the function documentation.   



### Find the optimum by dynamic programming

Bellman's algorithm to compute the optimal solution for all possible trajectories. 



```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
```





Plot the policy function:



```r
qplot(x_grid, x_grid - x_grid[opt$D[,1]], xlab="stock size", ylab="escapement")
```

![plot of chunk policyfn_plot](http://farm9.staticflickr.com/8023/7184657097_3fcb8d3470_o.png) 


and the value function:



```r
qplot(x_grid, opt$V, xlab="stock size", ylab="value")
```

![plot of chunk valuefn_plot](http://farm9.staticflickr.com/8154/7184657301_dc79a38b99_o.png) 







### Large growth error



```r
sigma_g <- 0.25    # Noise in population growth
sigma_m <- 0.01     # noise in stock assessment measurement
sigma_i <- 0.01     # noise in implementation of the quota
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
```






```r
SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps=19999)
```



```
Library ggplot2 loaded.
```



```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
```





Plot the policy function:



```r
qplot(x_grid, x_grid - x_grid[opt$D[,1]], xlab="stock size", ylab="escapement")
```

![plot of chunk policyfn_plot2](http://farm8.staticflickr.com/7071/7184723753_5e94866431_o.png) 


and the value function:



```r
qplot(x_grid, opt$V, xlab="stock size", ylab="value")
```

![plot of chunk valuefn_plot2](http://farm8.staticflickr.com/7215/7184723971_85608f5899_o.png) 


