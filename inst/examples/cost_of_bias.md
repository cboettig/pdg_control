






# Calculating the cost of bias  
 * author Carl Boettiger, <cboettig@gmail.com>
 * license: CC0


 * knitr-formatted [source code](https://github.com/cboettig/pdg_control/blob/master/inst/examples/cost_of_bias.Rmd)
 * [Cached data](http://two.ucdavis.edu/cboettig/data/cost_of_bias/)

Implements a numerical version of the SDP described in (Sethi _et. al._ 2005).  Then compute the optimal solution under different forms of uncertainty and compare the results.  





## Model setup 

We will assume a Beverton-Holt state equation / population dynamics function, <span> \( f(x,h) = \frac{A x}{1 + B x} \)</span>



```r
f <- BevHolt
pars <- c(1.5, 0.05)
K <- (pars[1] - 1)/pars[2]
```



with parameters A = `1.5` and B = `0.05`.  The positive stationary root of the model is given by <span>\( \frac{A-1}{B} \)</span>, `10`.   



```r
p <- 1
c0 <- 0.01
c1 <- 0
profit <- profit_harvest(price=p, c0 = c0, c1 = c1) 
```




We also assume a profit function of the form <span>\( \Pi = p h - \left( c_0  + c_1 \frac{h}{x} \right) \frac{h}{x} \)</span>, conditioned on <span>\( h > x \)</span> and <span>\(x > 0 \)</span>, with price p = `1`, c0 = `0.01`, and c1 = `0`.  




```r
x_grid <- seq(0, 2 * K, length = 100)  
h_grid <- x_grid  
```




and solve the problem on a discrete grid of `100` for stock size and range `0`, `20`.  We use the same set of gridpoints for the possible harvest levels. 


## Scenarios 

We calculate the stochastic transition matrix for the probability of going from any state \\(x_t \\) to any other state \\(x_{t+1}\\) the following year, for each possible choice of harvest \\( h_t \\).  This provides a look-up table for the dynamic programming calculations.

### No Uncertainty 

The first scenario considers the completely deterministic case.  



```r
sigma_g <- 0.0001    # Noise in population growth
z_g <- function() 1 
z_m <- function() 1 
z_i <- function() 1 
```






```r
deterministic_SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
```






```r
det_opt <- find_dp_optim(deterministic_SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
```




We simulate 100 replicates of this system.  We will used a fixed seed so that we can compare these replicates to simulations under different conditions.  (Of course the seed is irrelevant at this stage since this is actually deterministic).  



```r
set.seed(42)
sims_one <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0=K, det_opt$D, z_g, z_m, z_i, profit)
})
```



### Growth uncertainty 



```r
sigma_g <- 0.15    # Noise in population growth
z_g <- function() rlnorm(1,  0, sigma_g) # mean 1
z_m <- function() 1 
z_i <- function() 1 
```




The next scenario introduces growth uncertainty into the model, <span> \( x_{t+1} = z_g f(x_t) \) </span>, where `z_g` is lognormal with logmean 0 and logsd of `0.15`.  



```r
SDP_Mat <- determine_SDP_matrix(f, pars, x_grid, h_grid, sigma_g )
```






```r
opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
```





As before, we simulate 100 replicates using the same random number sequence, now under this case where the noise in growth is intrinsic and is being accounted for by the management.  



```r
set.seed(42)
sims_two <- lapply(1:100, function(i){
  ForwardSimulate(f, pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i, profit)
})
```




### Growth uncertainty & bias  




```r
est_pars <- pars
par_var <- .2
```




This time we consider the same optimization under uncertainty as before, but the simulations introduce bias through a random estimate of the growth rate parameter A, drawn from a normal with mean equal to the true value `1.5` and variance `0.2`.   Since A is a constant multiplier in the growth dynamics, this is equivalent to a random estimate of mean of the noise process `z_g`.  Our estimate of the parameter is drawn from the distribution and then held fixed for that replicate.  Each replicate draws its own value, so average parameter estimate across the replicates should be close to the true value.  _Isn't this equivalent to the standard parameter uncertainty?_



```r
bias <- rnorm(100, pars[1], par_var)
set.seed(42)
sims_three <- lapply(1:100, function(i){
  est_pars[1] <- bias[i]
  ForwardSimulate(f, est_pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i, profit)
})
```




For the record, the biases by index number are:



```r
bias
```



```
  [1] 1.766 1.380 1.511 1.394 1.484 1.532 1.480 1.429 1.406 1.507 1.460
 [12] 1.391 1.575 1.581 1.534 1.280 1.580 1.585 1.443 1.611 1.753 1.566
 [23] 1.653 1.757 1.832 1.404 1.418 1.434 1.720 1.616 1.692 1.560 1.442
 [34] 1.209 1.159 1.902 1.501 1.612 1.639 1.789 1.627 1.505 1.649 1.569
 [45] 1.374 1.692 1.533 1.493 1.316 1.730 1.422 1.480 1.494 1.489 1.611
 [56] 1.680 1.607 1.206 1.018 1.427 1.464 1.470 1.221 1.364 1.136 2.028
 [67] 1.609 1.499 1.122 1.447 1.664 1.656 1.591 1.782 1.138 1.280 1.496
 [78] 1.530 1.394 1.740 1.091 1.297 1.821 1.678 1.196 1.665 1.573 1.234
 [89] 1.572 1.679 1.301 1.604 1.253 1.606 1.303 1.425 1.292 1.640 1.471
[100] 1.600
```




## Known bias

The correct baseline comparison against the above scenario, Mike points out to me, still includes the bias, but solves for the optimal solution of that random draw.  This is much slower since it requires solving a new optimum each time, but parallelizes well.   




```r
#require(snowfall)
#sfInit(parallel=T, cpu=16)
#sfLibrary(pdgControl)
#sfExportAll()
set.seed(42)
sims_four <- lapply(1:100, function(i){
  est_pars[1] <- bias[i] 
  SDP_Mat <- determine_SDP_matrix(f, est_pars, x_grid, h_grid, sigma_g )
  opt <- find_dp_optim(SDP_Mat, x_grid, h_grid, OptTime=25, xT=0, 
                     profit, delta=0.05, reward=0)
  ForwardSimulate(f, est_pars, x_grid, h_grid, x0=K, opt$D, z_g, z_m, z_i, profit)
})
```






## Summarize and plot the results                                                   




```r
sims <- list(known = sims_one, Growth = sims_two, RandomBias = sims_three, KnownBias = sims_four)

dat <- melt(sims, id=names(sims_one[[1]]))  
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

![plot of chunk onerep](http://farm8.staticflickr.com/7079/7369822478_d53d7db9c6_o.png) 



This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again.



```r
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) 
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm8.staticflickr.com/7098/7369822864_441df2d76b_o.png) 





```r
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, profit))  + facet_wrap(~uncertainty)
```

![The profits made in each time interval of a single replicate, by scenario](http://farm8.staticflickr.com/7080/7369823126_07205d9e7e_o.png) 





```r
profits <-dt[ , sum(profit), by=c("reps", "uncertainty")] 
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm6.staticflickr.com/5113/7369823486_c489c1bdd5_o.png) 


Summary statistics 



```r
profits[, mean(V1), by=uncertainty]
```



```
     uncertainty    V1
[1,]       known 32.68
[2,]      Growth 34.34
[3,]  RandomBias 35.47
[4,]   KnownBias 36.98
```



```r
profits[, sd(V1), by=uncertainty]
```



```
     uncertainty     V1
[1,]       known  0.000
[2,]      Growth  3.963
[3,]  RandomBias 15.403
[4,]   KnownBias 16.346
```







```r
save(list="dt", file="bias.rda")
```




Direct comparisons: 



```r
setkey(dt, uncertainty)
cost_of_bias = dt[uncertainty=="KnownBias", sum(profit), by=reps]$V1 - dt[uncertainty=="RandomBias", sum(profit), by=reps]$V1
cost_of_bias
```



```
  [1]  4.805109  1.440203  0.348876  0.537600  0.108679 -0.030362  0.200340
  [8]  0.264686  0.310560  0.000000 -0.003407  0.134247  0.245258  0.242225
 [15]  0.001675  1.805395  0.625988  0.022581  0.311775  0.814401  2.470442
 [22]  0.848399  1.321014  2.034869  5.069526 -0.005955  0.799331 -0.005131
 [29]  0.216308 -0.125269  2.434977 -0.131973  0.400451  3.928943  3.021026
 [36]  7.311024  0.000000  1.557760  1.473729  4.662718  1.781252  0.000000
 [43]  0.992867 -0.203316  1.204456  1.022335  0.176032  0.000000  1.118997
 [50]  2.387117  0.157296 -0.001755  0.000000  0.000000  0.953200  1.116970
 [57]  0.208433  4.279790  3.942975  0.383384  0.028547  0.111130  3.012619
 [64]  0.678308  3.364385 19.830527  1.105768  0.000000  3.762469 -0.205451
 [71]  3.300783  1.577274 -0.675260  4.864242  3.611919  1.381712  0.000000
 [78]  0.405626  0.348840  2.235865  3.563851  2.994397  5.828749  0.262359
 [85]  2.191019  1.614311 -0.197514  2.931595 -0.372816  3.644167  1.783823
 [92] -0.161279  2.616421 -0.060540  1.433780 -0.011975  2.655753  2.166332
 [99] -0.001685  0.106752
```



```r
mean(cost_of_bias)
```



```
[1] 1.507
```



```r
sd(cost_of_bias)
```



```
[1] 2.46
```




# References

Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
"Fishery Management Under Multiple Uncertainty." _Journal of
Environmental Economics And Management_, *50*. ISSN 00950696,
<URL: http://dx.doi.org/10.1016/j.jeem.2004.11.005>.


