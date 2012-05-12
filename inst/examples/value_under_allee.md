




# Calculating the value of information
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0

 Implements a numerical version of the SDP described in 

```

Error in bibentry1(bibtype = "Article", textVersion = NULL, header = NULL,  : 
  A bibentry of bibtype 'Article' has to correctly specify the field(s): author

```

.
 Compute the optimal solution under different forms of uncertainty.  The true uncertainty is much smaller than that assumed, just a small uncertainty in growth.  We compare the case where the uncertainty is equal to this tiny growth uncertainty to policies derived under much greater uncertainty assumptions, and see how the resulting profit varies.   





Define noise parameters 


Chose the state equation / population dynamics function



```r
f <- Myer_harvest
pars <- c(1, 2, 6) 
p <- pars # shorthand 
K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 -0.5 * p[3] ) / 2
x0 <- K - sigma_g ^ 2 / 2 
```




Take a look at the model: 



```r
f
```



```
function (x, h, p) 
{
    sapply(x, function(x) max(0, p[1] * x^p[2]/(1 + x^p[2]/p[3]) - 
        h))
}
<environment: namespace:pdgControl>
```






We use a harvest-based profit function with default parameters



```r
profit <- profit_harvest(price=1, c0 = 0.01) 
```




The `profit_harvest` function has the form \\( \Pi = h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x} \\), conditioned on \\( h > x \\) and \\(x > 0 \\). 

Now we must set up the discrete grids for stock size and havest levels (which will use same resolution as for stock), in order to calculate the SDP solution.   Here we set the gridsize to 100.  



```r
x_grid <- seq(0, 2 * K, length = 100)  
h_grid <- x_grid  
```





# Scenarios: 

We calculate the stochastic transition matrix for the probability of going from any state \\(x_t \\) to any other state \\(x_{t+1}\\) the following year, for each possible choice of harvest \\( h_t \\).  This provides a look-up table for the dynamic programming calculations.  


In the Sethi case, computing the distribution over multiple sources of noise is actually quite difficult.  Simulation turns out to be more efficient than numerically integrating over each distribution.  


## No Uncertainty 




```r
sigma_g <- 0.01    # Noise in population growth
sigma_m <- 0.0     # noise in stock assessment measurement
sigma_i <- 0.0     # noise in implementation of the quota
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
```



Find the transition matrix 



```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
```




Find the optimum solution



```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
```




Simulate 



```r
sims_known <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i, profit)
})
```




## Growth uncertainty 




```r
sigma_g <- 0.5    # Noise in population growth
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
```



Find the transition matrix 



```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
```




Find the optimum solution



```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
```




Simulate.  The real world doesn't have as much uncertainty as we do in our knowledge of the state.  



```r
sigma_g <- 0.05    # Noise in population growth
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
sims_g <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0=K, opt$D, z_g, profit=profit)
})
```





## Growth & stock measurement uncertainty 




```r
sigma_g <- 0.5    # Noise in population growth
sigma_m <- 0.5     # noise in stock assessment measurement
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() 1
```




Find the transition matrix.  Use the simulation method to account for the extra uncertainties 



```r
require(snowfall) 
sfInit(parallel=TRUE, cpu=16)
```



```
R Version:  R version 2.14.1 (2011-12-22) 

```



```r
SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps=999)
```



```
Library ggplot2 loaded.
```




Find the optimum solution



```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
```




Simulate 



```r
sigma_g <- 0.05    # Noise in population growth
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() 1
z_i <- function() 1
sims_gm <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i, profit)
})
```







## Growth, stock measurement & implementation uncertainty 




```r
sigma_g <- 0.5    # Noise in population growth
sigma_m <- 0.5     # noise in stock assessment measurement
sigma_i <- 0.5     # noise in implementation of the quota
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() rlnorm(1,  0, sigma_m) # mean 1
z_i <- function() rlnorm(1,  0, sigma_i) # mean 1
```




Find the transition matrix.  Use the simulation method to account for the extra uncertainties 



```r
SDP_Mat <- SDP_by_simulation(f, pars, x_grid, h_grid, z_g, z_m, z_i, reps=999)
```



```
Library ggplot2 loaded.
```




Find the optimum solution



```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
```




Simulate 



```r
sigma_g <- 0.05    # Noise in population growth
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() 1
z_i <- function() 1
sims_gmi <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i, profit)
})
```






## Summarize and plot the results                                                   

R makes it easy to work with this big replicate data set.  We make data tidy (melt), fast (data.tables), and nicely labeled.



```r

sims <- list(known = sims_known, growth = sims_g, growth_stock = sims_gm, growth_stock_harvest = sims_gmi)

dat <- melt(sims, id=names(sims_known[[1]]))  
dt <- data.table(dat)
setnames(dt, c("L2", "L1"), c("reps", "uncertainty")) # names are nice
```




### Plots 

Let's begin by looking at the dynamics of a single replicate. The line shows Reed's S, the level above which the stock should be harvested (where catch should be the difference between stock and S).  To confirm that this policy is being followed, note that harvesting only occurs when the stock is above this line, and harvest is proportional to the amount by which it is above.  Change the replicate `reps==` to see the results from a different replicate.  



```r
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, fishstock)) +
  geom_abline(intercept=opt$S, slope = 0) +
  geom_line(aes(time, harvest), col="darkgreen") + 
  facet_wrap(~uncertainty) 
```

![plot of chunk onerep](http://farm6.staticflickr.com/5315/7179301258_23234a0a86_o.png) 



This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again.



```r
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) 
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.2) + facet_wrap(~uncertainty)
```

![plot of chunk all](http://farm8.staticflickr.com/7096/7179301716_89b6b82f12_o.png) 


We can also look at the harvest dynamics:



```r
p1 + geom_line(aes(time, harvest, group = reps), alpha = 0.1, col="darkgreen") + facet_wrap(~uncertainty)
```

![plot of chunk harvestplot](http://farm6.staticflickr.com/5323/7179302154_18a509be99_o.png) 


This strategy is supposed to be a constant-escapement strategy. We can visualize the escapement: 



```r
p1 + geom_line(aes(time, escapement, group = reps), alpha = 0.1, col="darkgrey") + facet_wrap(~uncertainty)
```

![plot of chunk escapement](http://farm6.staticflickr.com/5152/7179302548_bb0f02526a_o.png) 





```r
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, profit))  + facet_wrap(~uncertainty)
```

![plot of chunk unnamed-chunk-16](http://farm6.staticflickr.com/5076/7179303004_1195348c77_o.png) 




```r
profits <-dt[ , sum(profit), by=c("reps", "uncertainty")] 
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![plot of chunk unnamed-chunk-17](http://farm8.staticflickr.com/7099/7179303326_df650d1727_o.png) 


Summary stats



```r
profits[, mean(V1), by=uncertainty]
```



```
              uncertainty    V1
[1,]                known 18.07
[2,]               growth 17.87
[3,]         growth_stock 14.93
[4,] growth_stock_harvest 16.92
```



```r
profits[, sd(V1), by=uncertainty]
```



```
              uncertainty     V1
[1,]                known 0.1282
[2,]               growth 0.7367
[3,]         growth_stock 0.5289
[4,] growth_stock_harvest 0.6541
```






# References




