






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
sigma_g <- 0.0    # Noise in population growth
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
  ForwardSimulate(f, est_pars, x_grid, h_grid, x0=K, det_opt$D, z_g, z_m, z_i, profit)
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

![plot of chunk onerep](http://farm8.staticflickr.com/7088/7369448908_1e2c1752ae_o.png) 



This plot summarizes the stock dynamics by visualizing the replicates. Reed's S shown again.



```r
p1 <- ggplot(dt) + geom_abline(intercept=opt$S, slope = 0) 
p1 + geom_line(aes(time, fishstock, group = reps), alpha = 0.1) + facet_wrap(~uncertainty)
```

![the induced dynamics in the stock size over time, for all replicates, by scenario](http://farm9.staticflickr.com/8151/7184215351_ca5ae88525_o.png) 





```r
ggplot(subset(dt,reps==1)) +
  geom_line(aes(time, profit))  + facet_wrap(~uncertainty)
```

![The profits made in each time interval of a single replicate, by scenario](http://farm6.staticflickr.com/5192/7184215609_d9f30778e3_o.png) 





```r
profits <-dt[ , sum(profit), by=c("reps", "uncertainty")] 
ggplot(profits) + geom_histogram(aes(V1)) + facet_wrap(~uncertainty)
```

![the distribution of profits by scenario](http://farm8.staticflickr.com/7227/7369450174_1ee9f68bc8_o.png) 


Summary statistics 



```r
profits[, mean(V1), by=uncertainty]
```



```
     uncertainty    V1
[1,]       known 10.03
[2,]      Growth 34.34
[3,]  RandomBias 35.47
[4,]   KnownBias 10.03
```



```r
profits[, sd(V1), by=uncertainty]
```



```
     uncertainty       V1
[1,]       known  0.00000
[2,]      Growth  3.96332
[3,]  RandomBias 15.40344
[4,]   KnownBias  0.05074
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
  [1] -50.28123  -9.62648 -29.31469 -17.39949 -21.86984 -20.82687 -19.43285
  [8] -17.54376 -16.80563 -24.01193 -23.07597 -19.23580 -30.31750 -30.49961
 [15] -25.09872  -7.76515 -33.52429 -27.95587 -17.82549 -36.02269 -43.17873
 [22] -35.92103 -37.45331 -41.79552 -50.51585 -20.04248 -12.98206 -21.03775
 [29] -35.34005 -29.04389 -42.35734 -24.44153 -17.63498  -0.06859   0.91319
 [36] -56.38475 -29.91497 -40.39306 -38.68930 -48.83867 -39.52322 -34.13575
 [43] -36.10616 -25.66183 -10.94967 -37.14417 -30.90066 -23.79696 -10.95664
 [50] -42.58770 -20.45485 -22.30158 -23.06129 -24.07780 -35.96038 -37.06648
 [57] -29.87760   1.17896   4.13925 -18.02250 -19.97941 -20.03482  -3.08429
 [64] -14.90505   3.17125 -79.31692 -37.09696 -24.82952   0.30099 -25.67112
 [71] -47.16303 -37.96391 -23.11786 -50.06213   1.99268  -8.98098 -27.43131
 [78] -32.35128 -17.78517 -41.80647   3.93866  -5.15093 -52.56687 -32.64158
 [85]  -5.44303 -38.80474 -25.52556  -1.69115 -25.24099 -46.67603  -9.93402
 [92] -28.24107  -6.04543 -27.71433 -10.09559 -22.90526  -6.51199 -40.96808
 [99] -24.87223 -29.82843
```



```r
mean(cost_of_bias)
```



```
[1] -25.44
```



```r
sd(cost_of_bias)
```



```
[1] 15.37
```




# References

Sethi G, Costello C, Fisher A, Hanemann M and Karp L (2005).
"Fishery Management Under Multiple Uncertainty." _Journal of
Environmental Economics And Management_, *50*. ISSN 00950696,
<URL: http://dx.doi.org/10.1016/j.jeem.2004.11.005>.


